t1 = time()
import MUMPS
import MPI
import SparseArrays
include("../core/SystemReader.jl")
include("../core/SystemExporter.jl")
t2 = time()
println("Time taken to import modules: $(t2 - t1) seconds")

""" Solve system of equations using MUMPS.

This function is run on multiple CPUs. Only the CPU with rank 0 reads
the system of equations and saves the resulting solution vector to file.

# Arguments
- `soe_dir_path`: Path to system of equation directory containing the system files
- `nthreads`: Number of threads per core for OpenMP
- `analysis_method_itype`: Analysis step method (1 = serial, 2 = parallel)
- `parallel_ordering_method_itype`: Parallel ordering method (1 = PT-SCOTCH, 2 = ParMETIS)
- `verbose_output_itype`: Verbose output flag (0 = silent, 1 = verbose)
- `memory_relax_perc`: Memory relaxation percentage
"""
function solve_system(
    soe_dir_path::String,
    nthreads::Int,
    analysis_method_itype::Int,
    parallel_ordering_method_itype::Int,
    verbose_output_itype::Int,
    memory_relax_perc::Int
)::Nothing
    MPI.Init()
    comm = MPI.COMM_WORLD
    root = 0
    mumps = MUMPS.Mumps{Float64}(MUMPS.mumps_unsymmetric, MUMPS.default_icntl, MUMPS.default_cntl64)
    if MPI.Comm_rank(comm) == root
        # Set analysis step parameters
        MUMPS.set_icntl!(mumps, 28, analysis_method_itype)  # Analysis step method
        MUMPS.set_icntl!(mumps, 29, parallel_ordering_method_itype)  # Parallel ordering method
        MUMPS.set_icntl!(mumps, 14, memory_relax_perc)  # Memory relaxation percentage
        (N, Li, Lj, Lv, rhs) = SystemReader.read_system_of_equations_from_file(soe_dir_path)
        A = SparseArrays.sparse(Li, Lj, Lv, N, N)
        MUMPS.associate_matrix!(mumps, A)
        MUMPS.associate_rhs!(mumps, rhs)
        if verbose_output_itype == 1
            print_icntl_parameters(mumps)
        end
    end
    MUMPS.factorize!(mumps)
    MUMPS.solve!(mumps)
    MPI.Barrier(comm)
    if MPI.Comm_rank(comm) == root
        x = vec(MUMPS.get_solution(mumps))
        SystemExporter.export_solution_vector(soe_dir_path, x)
    end
    MUMPS.finalize(mumps)
    MPI.Finalize()
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

if abspath(PROGRAM_FILE) == @__FILE__
    soe_dir_path = ARGS[1]
    src_dir_path = ARGS[2]
    nthreads = parse(Int, ARGS[3])
    analysis_method = ARGS[4]
    parallel_ordering_method = ARGS[5]
    verbose_output_itype = parse(Int, ARGS[6])
    memory_relax_perc = parse(Int, ARGS[7])
    (
        analysis_method_itype, 
        parallel_ordering_method_itype
    ) = get_integer_option_values(
        analysis_method,
        parallel_ordering_method
    )
    solve_system(
        soe_dir_path,
        nthreads,
        analysis_method_itype,
        parallel_ordering_method_itype,
        verbose_output_itype,
        memory_relax_perc
    )
end 