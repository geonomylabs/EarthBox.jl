module ParallelSolver

include("core/NamesManager.jl")
include("core/SystemReader.jl")
include("core/SystemExporter.jl")
include("core/MumpsChecks.jl")
include("core/MumpsLoop.jl")
include("run_mumps/RunMumpsSolverExternal.jl")
include("run_mumps/RunMumpsSolverLoop.jl")

import EarthBox.PrintFuncs: print_info
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState
import EarthBox.ConfigurationManager.SolverConfig: modify_mumps_parameters!
import EarthBox.ConfigurationManager.SolverConfig: print_mumps_solver_info
import EarthBox.ConfigurationManager.SolverConfig: reset_mumps_parameters_to_original!
import .SystemReader: read_solution_vector_file
import .SystemExporter: export_system_of_equations_for_mumps
import .RunMumpsSolverExternal
import .RunMumpsSolverLoop


"""
    parallel_direct_solver(N::Int64, Li::Vector{Int64}, Lj::Vector{Int64}, Lv::Vector{Float64}, 
        R::Vector{Float64}, solver_config::SolverConfigState)::Vector{Float64})::Vector{Float64}

Solve system of equations using MUMPS solver.

A variety of options are available for how MUMPS is used:
- Option 1: A script with a loop around the mumps solver is spawned as a child process and communicates with 
   the parent process via file IO. This is the default behavior. This option is activated if `use_internal_mumps`
   is set to true and `pass_large_arrays_via_mpi` is set to false in the solver configuration.
- Option 2: A script with a loop around the mumps solver is spawned as a child process and communicates with 
    the parent process via MPI Send and Receive. This approach does not work on all systems and is therefore 
    tuned off by default. This option is activated if `use_internal_mumps` and `pass_large_arrays_via_mpi` 
    are both set to true in the solver configuration.
- Option 3: The mumps solver is run via the command line using the `run()` function and communicates with 
    the parent process via file IO. This option is activated if `use_internal_mumps` is set to false in the 
    solver configuration.

# Arguments
- `N::Int64`: Number of unknowns
- `Li::Vector{Int64}`: Large matrix 1-based row indices for non-zero matrix elements
- `Lj::Vector{Int64}`: Large matrix 1-based column indices for non-zero matrix elements
- `Lv::Vector{Float64}`: Non-zero matrix values of discretized system of equations
- `R::Vector{Float64}`: Right-hand-side vector of discretized equations
- `solver_config::SolverConfigState`: Solver configuration parameters

# Returns
- `Vector{Float64}`: Solution vector
"""
function parallel_direct_solver(
    N::Int64,
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    R::Vector{Float64},
    solver_config::SolverConfigState
)::Vector{Float64}

    output_dir         = solver_config.output_dir
    src_dir            = solver_config.src_dir
    nproc_orig         = solver_config.nproc
    ordering_orig      = solver_config.parallel_ordering_method
    memory_relax_orig  = solver_config.memory_relax_perc
    use_internal_mumps = solver_config.use_internal_mumps
    mpi_comm           = solver_config.mpi_comm
    pass_large_arrays_via_mpi = solver_config.pass_large_arrays_via_mpi

    if !pass_large_arrays_via_mpi
        soe_dir_path = export_system_of_equations_for_mumps(output_dir, N, Li, Lj, Lv, R)
    end
    
    nmumps_max = 5
    nmumps = 0
    mumps_solver_status = "Failure"
    S = zeros(Float64, 2)
    while mumps_solver_status == "Failure" && nmumps < nmumps_max
        if nmumps == 0
            nproc_orig = solver_config.nproc
            ordering_orig = solver_config.parallel_ordering_method
            memory_relax_orig = solver_config.memory_relax_perc
        end
        modify_mumps_parameters!(solver_config, nmumps)
        print_mumps_solver_info(solver_config, nmumps)
        if mpi_comm !== nothing && use_internal_mumps
            matrix_info = [N, length(Li)]
            if pass_large_arrays_via_mpi
                S, mumps_solver_status = RunMumpsSolverLoop.mpirun_mumps_mpi_comm(
                    solver_config, Li, Lj, Lv, R, matrix_info)
            else
                termination_flag = 0
                mumps_solver_status_io_comm = RunMumpsSolverLoop.mpirun_mumps_io_comm(
                    solver_config, termination_flag)
                print_info("Solver status from IO communication: $(mumps_solver_status_io_comm)", level=2)
                S, mumps_solver_status = read_solution_vector_file(soe_dir_path)
            end
        else
            RunMumpsSolverExternal.mpirun_mumps(output_dir, src_dir, solver_config)
            S, mumps_solver_status = read_solution_vector_file(soe_dir_path)
        end
        print_solver_status(mumps_solver_status)
        reset_mumps_parameters_to_original!(
            solver_config, nproc_orig, ordering_orig, memory_relax_orig)
        nmumps += 1
    end
    if mumps_solver_status == "Failure"
        error("MUMPS solver failed after multiple attempts.")
    end
    return S
end

function print_solver_status(mumps_solver_status::String)
    description = "MUMPS solver status: "
    msg = "$(description) $(mumps_solver_status)"
    print_info(msg, level=2)
end

end # module ParallelSolver 