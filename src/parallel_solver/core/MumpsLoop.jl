module MumpsLoop

t1 = time()
import MUMPS
import MPI
import SparseArrays
import Base.catch_backtrace
import LinearAlgebra: norm
include("./SystemReader.jl")
include("./SystemExporter.jl")
include("./NamesManager.jl")
t2 = time()
println("Time taken to import modules: $(t2 - t1) seconds")

function print_info(msg::String; level::Int=1)
    if level == 1
        println("<*> $msg")
    elseif level == 2
        println("   <*> $msg")
    elseif level == 3
        println("      <*> $msg")
    end
end

function print_error(msg::String; level::Int=1)
    if level == 1
        println("<*> !!!ERROR!!! $msg")
    elseif level == 2
        println("   <*> !!!ERROR!!! $msg")
    elseif level == 3
        println("      <*> !!!ERROR!!! $msg")
    end
end

""" Solve system of equations using MUMPS.

This function is run on multiple CPUs. Only the CPU with rank 0 reads
the system of equations and saves the resulting solution vector to file.

# Arguments
- `analysis_method_itype`: Analysis step method (1 = serial, 2 = parallel)
- `parallel_ordering_method_itype`: Parallel ordering method (1 = PT-SCOTCH, 2 = ParMETIS)
- `verbose_output_itype`: Verbose output flag (0 = silent, 1 = verbose)
- `memory_relax_perc`: Memory relaxation percentage
"""
function solve_system_io_comm(
    analysis_method_itype::Int,
    parallel_ordering_method_itype::Int,
    verbose_output_itype::Int,
    memory_relax_perc::Int,
    soe_dir_path::String
)::Nothing
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    root = 0
    loop_check_time = 5.0
    if rank == root
        print_info("Starting MUMPS solver loop with large arrays passed via IO on rank $rank", level=2)
        print_info("SOE directory path: $soe_dir_path", level=2)
    end
    if !isdir(soe_dir_path)
        mkpath(soe_dir_path)
        print_info("Created directory: $soe_dir_path", level=2)
    end
    try
        start_print_time = time()
        last_print_time = start_print_time
        while true
            total_time = time() - start_print_time
            delta_time = time() - last_print_time
            ready_for_solver = false
            if rank == root
                ready_for_solver = has_a_ready_to_solve_file(soe_dir_path)
                MPI.Bcast!(Ref(ready_for_solver), root, comm)
            else
                ready_for_solver = Ref{Bool}()
                MPI.Bcast!(ready_for_solver, root, comm)
                ready_for_solver = ready_for_solver[]
            end
            if rank == root && delta_time >= loop_check_time
                print_info("MUMPS solver loop with file IO communication is active on rank $rank for $total_time seconds with ready_for_solver: $ready_for_solver", level=2)
                last_print_time = time()
            end
            if !ready_for_solver
                sleep(0.01)
            else
                mumps = MUMPS.Mumps{Float64}(
                    MUMPS.mumps_unsymmetric, MUMPS.default_icntl, MUMPS.default_cntl64)
                # Set analysis step parameters on all processes
                if rank == root
                    MUMPS.set_icntl!(mumps, 28, analysis_method_itype)
                    MUMPS.set_icntl!(mumps, 29, parallel_ordering_method_itype)
                    MUMPS.set_icntl!(mumps, 14, memory_relax_perc)
                    if verbose_output_itype == 0
                        MUMPS.set_icntl!(mumps, 1, 6)  # Error messages to stderr
                        MUMPS.set_icntl!(mumps, 2, 0)  # Diagnostics suppressed
                        MUMPS.set_icntl!(mumps, 3, 0)  # Global info suppressed
                        MUMPS.set_icntl!(mumps, 4, 0)  # Print level suppressed
                    end

                    # Remove the ready to solve file since it has been read
                    remove_ready_to_solve_file(soe_dir_path)
                    # Read the termination info from the termination file
                    termination_flag = manage_termination_info(soe_dir_path)
                    # Check for termination signal (negative N means stop)
                    if termination_flag < 0
                        print_info("Termination signal received. Stopping loop.", level=2)
                        break
                    end
                    
                    (N, Li, Lj, Lv, rhs) = SystemReader.read_system_of_equations_from_file(soe_dir_path)

                    A = SparseArrays.sparse(Li, Lj, Lv, N, N)
                    MUMPS.associate_matrix!(mumps, A)
                    MUMPS.associate_rhs!(mumps, rhs)

                    if verbose_output_itype == 1
                        print_icntl_parameters(mumps)
                    end
                end
                
                MUMPS.factorize!(mumps)
                if rank == root
                    found_error_factorization = check_mumps_errors(mumps, "factorization", rank)
                end
                
                MUMPS.solve!(mumps)
                if rank == root
                    found_error_solve = check_mumps_errors(mumps, "solve", rank)
                end

                MPI.Barrier(comm)
                if rank == root
                    x = vec(MUMPS.get_solution(mumps))
                    solution_flag = 0
                    if !found_error_solve
                        SystemExporter.export_solution_vector(soe_dir_path, x)
                        solution_flag = 1
                    else
                        solution_flag = 0
                    end
                    # Communicate with the parent process that the solve processes is done by
                    # creating a solution flag file that contains the solution flag indicating
                    # success or failure of the solve process
                    create_solution_flag_file(soe_dir_path, solution_flag)
                end
                MUMPS.finalize(mumps)
            end
        end
    catch e
        if rank == root
            error_msg = "Mumps solve_system error: $e"
            print_error("rank $rank : $error_msg", level=2)
            create_error_flag_file(soe_dir_path, 0)
        end
    finally
        MPI.Finalize()
        print_info("MUMPS solver loop finalized on rank $rank", level=2)
    end
    return nothing
end

function manage_termination_info(soe_dir_path::String)::Int64
    names = NamesManager.FileAndDirNames()
    file_name = names.termination_info_file_name
    file_path = joinpath(soe_dir_path, file_name)
    termination_flag = open(file_path, "r") do f
        read(f, Int64)
    end
    rm(file_path)
    return termination_flag
end

function remove_ready_to_solve_file(soe_dir_path::String)::Nothing
    names = NamesManager.FileAndDirNames()
    file_name = names.ready_to_solve_file_name
    file_path = joinpath(soe_dir_path, file_name)
    if isfile(file_path)
        rm(file_path)
    end
    return nothing
end

function has_a_ready_to_solve_file(soe_dir_path::String)::Bool
    names = NamesManager.FileAndDirNames()
    file_name = names.ready_to_solve_file_name
    file_path = joinpath(soe_dir_path, file_name)
    return isfile(file_path)
end

function create_solution_flag_file(soe_dir_path::String, solution_flag::Int64)::Nothing
    names = NamesManager.FileAndDirNames()
    solution_file_name = names.solution_flag_file_name
    solution_file_path = joinpath(soe_dir_path, solution_file_name)
    open(solution_file_path, "w") do f
        write(f, solution_flag)
    end
    return nothing
end

function get_ready_to_solve_files(soe_dir_path::String)::Vector{String}
    flag_files = filter(f -> startswith(f, "ready_to_solve_"), readdir(soe_dir_path))
    return flag_files
end

function get_termination_info_files(soe_dir_path::String)::Vector{String}
    termination_info_files = filter(f -> startswith(f, "termination_info_"), readdir(soe_dir_path))
    return termination_info_files
end

function solve_system_mpi_comm(
    analysis_method_itype::Int,
    parallel_ordering_method_itype::Int,
    verbose_output_itype::Int,
    memory_relax_perc::Int,
    soe_dir_path::String
)::Nothing
    MPI.Init()
    parent = MPI.Comm_get_parent()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    root = 0
    loop_check_time = 5.0
   
    print_info("Starting MUMPS solver loop with large arrays passed via MPI on rank $rank", level=2)
    if rank == root
        print_info("SOE directory path: $soe_dir_path", level=2)
    end
    try
        start_print_time = time()
        last_print_time = start_print_time
        while true
            total_time = time() - start_print_time
            delta_time = time() - last_print_time
            if rank == root
                has_matrix_info, status = MPI.Iprobe(0, 0, parent)
                MPI.Bcast!(Ref(has_matrix_info), root, comm)
            else
                has_matrix_info = Ref{Bool}()
                MPI.Bcast!(has_matrix_info, root, comm)
                has_matrix_info = has_matrix_info[]
            end
            #if rank == root && delta_time >= loop_check_time
            if delta_time >= loop_check_time
                print_info("MUMPS solver loop with MPI array passing active on rank $rank for $total_time seconds with has_matrix_info: $has_matrix_info", level=2)
                last_print_time = time()
            end
            if has_matrix_info
                mumps = MUMPS.Mumps{Float64}(
                    MUMPS.mumps_unsymmetric, MUMPS.default_icntl, MUMPS.default_cntl64)
                # Set analysis step parameters on all processes
                if rank == root
                    MUMPS.set_icntl!(mumps, 28, analysis_method_itype)
                    MUMPS.set_icntl!(mumps, 29, parallel_ordering_method_itype)
                    MUMPS.set_icntl!(mumps, 14, memory_relax_perc)
                    if verbose_output_itype == 0
                        MUMPS.set_icntl!(mumps, 1, 6)   # Error messages to stderr
                        MUMPS.set_icntl!(mumps, 2, 0)  # Diagnostics suppressed
                        MUMPS.set_icntl!(mumps, 3, 0)  # Global info suppressed
                        MUMPS.set_icntl!(mumps, 4, 0)   # Print level suppressed
                    end
                    matrix_info = Vector{Int64}(undef, 2)

                    MPI.Recv!(matrix_info, 0, 0, parent)

                    N = matrix_info[1]
                    nnz = matrix_info[2]
                    # Check for termination signal (negative N means stop)
                    if N < 0
                        print_info("Termination signal received. Stopping loop.", level=2)
                        break
                    end
                    
                    Li = Vector{Int64}(undef, nnz)
                    Lj = Vector{Int64}(undef, nnz)
                    Lv = Vector{Float64}(undef, nnz)
                    rhs = Vector{Float64}(undef, N)

                    MPI.Recv!(Li, 0, 1, parent)
                    MPI.Recv!(Lj, 0, 2, parent)
                    MPI.Recv!(Lv, 0, 3, parent)
                    MPI.Recv!(rhs, 0, 4, parent)

                    A = SparseArrays.sparse(Li, Lj, Lv, N, N)
                    MUMPS.associate_matrix!(mumps, A)
                    MUMPS.associate_rhs!(mumps, rhs)

                    if verbose_output_itype == 1
                        print_icntl_parameters(mumps)
                    end
                end
                
                MUMPS.factorize!(mumps)
                if rank == root
                    found_error_factorization = check_mumps_errors(mumps, "factorization", rank)
                end
                
                MUMPS.solve!(mumps)
                if rank == root
                    found_error_solve = check_mumps_errors(mumps, "solve", rank)
                end

                MPI.Barrier(comm)
                if rank == root
                    x = vec(MUMPS.get_solution(mumps))
                    if found_error_solve == true
                        x[1] = NaN
                    end
                    print_info("Sending solution to parent...", level=2)
                    MPI.Send(x, 0, 5, parent)
                    print_info("Solution sent to parent...", level=2)
                end
                MUMPS.finalize(mumps)
            end
        end
    catch e
        if rank == root
            error_msg = "Mumps solve_system error: $e"
            print_error("rank $rank : $error_msg", level=2)
        end
    finally
        MPI.Finalize()
        print_info("MUMPS solver loop finalized on rank $rank", level=2)
    end
    return nothing
end

function manage_error(rank::Int, error_msg::String, parent::MPI.Comm)::Nothing
    print_error("rank $rank : $error_msg", level=2)
    error_flag = [1]
    print_info("Sending error flag to parent...", level=2)
    MPI.Send(error_flag, 0, 6, parent)
    print_info("Error flag sent to parent...", level=2)
    return nothing
end

""" Print the integer control parameters for the MUMPS solver.
"""
function print_icntl_parameters(mumps::MUMPS.Mumps)
    println("\nInteger Control (icntl) Parameters for MUMPS Solver")
    println(">> icntl 28 = ", mumps.icntl[28], ": 1 = serial analysis step; 2 = parallel")
    println(">> icntl 29 = ", mumps.icntl[29], ": 1 = PT-SCOTCH; 2 = ParMetis")
    println(">> icntl 14 = ", mumps.icntl[14], ": Percentage increase in the estimated working space")
    println(">> icntl  8 = ", mumps.icntl[8], ": Scaling strategy: 0 = off; 77 = automatic")
    println(">> icntl 23 = ", mumps.icntl[23], ": Maximum size of working memory allocated per core")
end

function get_integer_option_values(
    analysis_method::String,
    parallel_ordering_method::String
)::Tuple{Int, Int}
    analysis_method_itype = analysis_method == "PARALLEL" ? 2 : 1
    parallel_ordering_method_itype = parallel_ordering_method == "PTSCOTCH" ? 1 : 2
    return (analysis_method_itype, parallel_ordering_method_itype)
end

function manage_mumps_error(mumps::MUMPS.Mumps, phase::String)
    # info1 = mumps.info[1] # error with largest magnitude (same on all processes) (-1 = error, 0 = no error/warning, +1 = warning)
    # info2 = mumps.info[2] # error information
    infog1 = mumps.infog[1] # global error with largest magnitude (same on all processes) (-1 = error, 0 = no error/warning, +1 = warning)
    infog2 = mumps.infog[2] # error/warning information
    if infog1 < 0
        error_msg = "MUMPS error in $phase phase: INFOG(1)=$infog1, INFOG(2)=$infog2"
        throw(ErrorException(error_msg))
    end
end

function check_mumps_errors(mumps::MUMPS.Mumps, phase::String, rank::Int)::Bool
    info1 = mumps.info[1] # error with largest magnitude (same on all processes) (-1 = error, 0 = no error/warning, +1 = warning)
    info2 = mumps.info[2] # error information
    infog1 = mumps.infog[1] # global error with largest magnitude (same on all processes) (-1 = error, 0 = no error/warning, +1 = warning)
    infog2 = mumps.infog[2] # error/warning information
    found_error = false
    if infog1 < 0
        error_msg = "MUMPS error on rank $rank in $phase phase: INFOG(1)=$infog1, INFOG(2)=$infog2"
        println("   >> $error_msg")
        found_error = true
    end
    if info1 < 0
        error_msg = "MUMPS error on rank $rank in $phase phase: INFO(1)=$info1, INFO(2)=$info2"
        println("   >> $error_msg")
        found_error = true
    end
    return found_error
end

end # module