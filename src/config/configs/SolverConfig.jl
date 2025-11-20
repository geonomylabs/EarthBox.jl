module SolverConfig

import YAML
import MPI
import EarthBox.PrintFuncs: print_info

Base.@kwdef mutable struct InternalMumpsSolver
    use_mumps::Bool = false
    inter::Union{MPI.Comm, Nothing} = nothing
    is_running::Bool = false
    comm::Union{MPI.Comm, Nothing} = nothing
end

"""
    SolverConfigState

Mutable struct to store solver configuration parameters.

# Fields
- `use_mumps::Bool`: Use MUMPS solver
- `nthreads::Int`: Number of threads
- `nproc::Int`: Number of processors
- `analysis_method::String`: Analysis method (PARALLEL or SERIAL)
- `parallel_ordering_method::String`: Parallel ordering method (PTSCOTCH or ParMETIS)
- `memory_relax_perc::Int`: Memory relaxation percentage
- `verbose_output::Int`: Verbose output level
- `pymumps_timeout::Float64`: Timeout for MUMPS solver in seconds
- `output_dir::Union{String, Nothing}`: Output directory for system files
- `src_dir::Union{String, Nothing}`: Source directory for parallel solver package
"""
mutable struct SolverConfigState
    use_mumps::Bool
    nthreads::Int
    nproc::Int
    use_internal_mumps::Bool
    analysis_method::String
    parallel_ordering_method::String
    memory_relax_perc::Int
    verbose_output::Int
    pymumps_timeout::Float64
    output_dir::Union{String, Nothing}
    src_dir::Union{String, Nothing}
    mpi_comm::Union{MPI.Comm, Nothing}
    mpi_initialized::Bool
    mpi_rank::Int
    pass_large_arrays_via_mpi::Bool
    internal_mumps_solver::InternalMumpsSolver
end

function SolverConfigState(;
    use_mumps::Bool = false,
    nthreads::Int = 1,
    nproc::Int = 1,
    analysis_method::String = "SERIAL",
    parallel_ordering_method::String = "ParMETIS",
    memory_relax_perc::Int = 25,
    verbose_output::Int = 0,
    pymumps_timeout::Float64 = 3600.0,
    output_dir::Union{String, Nothing} = nothing,
    src_dir::Union{String, Nothing} = nothing,
    use_internal_mumps::Bool = false,
    mpi_comm::Union{MPI.Comm, Nothing} = nothing,
    mpi_initialized::Bool = false,
    mpi_rank::Int = 0,
    pass_large_arrays_via_mpi::Bool = false
)
    return SolverConfigState(
        use_mumps,
        nthreads,
        nproc,
        use_internal_mumps,
        analysis_method,
        parallel_ordering_method,
        memory_relax_perc,
        verbose_output,
        pymumps_timeout,
        output_dir,
        src_dir,
        mpi_comm,
        mpi_initialized,
        mpi_rank,
        pass_large_arrays_via_mpi,
        InternalMumpsSolver(use_mumps=use_mumps)
    )
end

function update_output_directory!(config::SolverConfigState, output_dir::String)
    config.output_dir = output_dir
end

function update_src_directory!(config::SolverConfigState, src_dir::String)
    config.src_dir = src_dir
end

function read_solver_config!(config::SolverConfigState, solver_options_file::String)
    solver_config = YAML.load_file(solver_options_file)
    check_solver_config!(solver_config)
    config.use_mumps = solver_config["use_mumps"]
    config.nthreads = solver_config["nthreads"]
    config.nproc = solver_config["nproc"]
    config.analysis_method = solver_config["analysis_method"]
    config.parallel_ordering_method = solver_config["parallel_ordering_method"]
    config.memory_relax_perc = solver_config["memory_relax_perc"]
    config.verbose_output = solver_config["verbose_output"]
end

function modify_mumps_parameters!(config::SolverConfigState, nmumps::Int)
    if nmumps == 1
        print_condition1_info()
    elseif nmumps == 2
        print_condition2_info()
        adjust_memory_relax!(config)
    elseif nmumps == 3
        set_alternative_parallel_ordering_method!(config)
        print_condition3_info(config.parallel_ordering_method)
    elseif nmumps == 4
        set_alternative_parallel_ordering_method!(config)
        update_nproc!(config)
        print_condition4_info(config.parallel_ordering_method, config.nproc)
    end
end

function adjust_memory_relax!(config::SolverConfigState, memory_relax_perc::Int = 35)
    config.memory_relax_perc = memory_relax_perc
end

function set_alternative_parallel_ordering_method!(config::SolverConfigState)
    config.parallel_ordering_method = config.parallel_ordering_method == "PTSCOTCH" ? "ParMETIS" : "PTSCOTCH"
end

function update_nproc!(config::SolverConfigState)
    config.nproc = max(1, config.nproc - 1)
end

function reset_mumps_parameters_to_original!(
    config::SolverConfigState,
    nproc_orig::Int,
    ordering_orig::String,
    memory_relax_orig::Int
)
    config.nproc = nproc_orig
    config.parallel_ordering_method = ordering_orig
    config.memory_relax_perc = memory_relax_orig
end

function print_mumps_solver_info(config::SolverConfigState, nmumps::Int)
    print_info(
        "MUMPS solver run $(nmumps): number of processors $(config.nproc): "
        *"analysis_method $(config.analysis_method): parallel_ordering_method $(config.parallel_ordering_method): "
        *"memory_relax_perc $(config.memory_relax_perc)",
        level=2
    )
end

function check_solver_config!(solver_config::Dict)
    if !haskey(solver_config, "nproc")
        solver_config["nproc"] = 2
    end
    if !haskey(solver_config, "nthreads")
        solver_config["nthreads"] = 1
    end
    if !haskey(solver_config, "analysis_method")
        solver_config["analysis_method"] = "PARALLEL"
    end
    if !haskey(solver_config, "parallel_ordering_method")
        solver_config["parallel_ordering_method"] = "ParMETIS" # "PTSCOTCH"
    end
    if !haskey(solver_config, "verbose_output")
        solver_config["verbose_output"] = 0
    end
    if !haskey(solver_config, "memory_relax_perc")
        solver_config["memory_relax_perc"] = 25
    end
end

function print_condition1_info()
    description = "Try again with the same parameters"
    print_info(description, level=3)
end

function print_condition2_info()
    description = "Try again with memory relaxation set to 35%"
    print_info(description, level=3)
end

function print_condition3_info(parallel_ordering_method::String)
    description = "Try again with parallel ordering method: $(parallel_ordering_method)"
    print_info(description, level=3)
end

function print_condition4_info(parallel_ordering_method::String, nproc::Int)
    description1 = "Try again with parallel ordering method:"
    description2 = "and one less processor, nprocs"
    msg = "$(description1) $(parallel_ordering_method) $(description2) $(nproc)"
    print_info(msg, level=3)
end

end # module SolverConfig 