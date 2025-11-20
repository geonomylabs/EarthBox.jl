module MPIManager

import MPI
import EarthBox.PrintFuncs: print_info
import MPIPreferences
# Set to true during code testing in the Julia REPL to avoid MPI finalization
# and allow rerunning the code without restarting Julia.
const DEBUG_MODE = false

function initialize_mpi!(use_mumps::Bool, use_internal_mumps::Bool)::Union{MPI.Comm, Nothing}
    if use_mumps && use_internal_mumps
        if !MPI.Initialized()
            MPI.Init()
            comm = MPI.COMM_WORLD
            rank = MPI.Comm_rank(comm)
            size = MPI.Comm_size(comm)
            if size > 1
                error(get_size_error_msg())
            end
            print_info("MPI initialized successfully on rank $rank of $size")
        else
            comm = MPI.COMM_WORLD
            rank = MPI.Comm_rank(comm)
            size = MPI.Comm_size(comm)
            if size > 1
                error(get_size_error_msg())
            end
            print_info("MPI already initialized on rank $rank of $size")
        end
    else
        comm = nothing
    end
    return comm
end

function get_size_error_msg()::String
    return (
        "Parallel EarthBox calculations must be run on a single parent processor. " 
        *"The parent processor spawns child processes to run the parallel MUMPS solver. "
        *"Use the `nprocs` argument to specify the number of child processes to spawn."
        )
end

function finalize_mpi!()::Nothing
    if MPI.Initialized()
        # Do not finalize MPI if in DEBUG_MODE
        if !DEBUG_MODE
            MPI.Finalize()
            print_info("MPI finalized successfully")
        end
    else
        print_info("MPI was not initialized")
    end
    return nothing
end

function is_mpi_initialized()::Bool
    return MPI.Initialized()
end

function get_mpi_comm()::Union{MPI.Comm, Nothing}
    if MPI.Initialized()
        return MPI.COMM_WORLD
    else
        return nothing
    end
end

function get_mpi_rank(comm::Union{MPI.Comm, Nothing})::Int
    if comm !== nothing
        return MPI.Comm_rank(comm)
    else
        return 0
    end
end

"""
    get_mpi_implementation()

Returns information about the MPI implementation being used.
Returns a NamedTuple with `binary` and `is_system_mpi` fields.
"""
function get_mpi_implementation()::NamedTuple
    try
        binary = MPIPreferences.binary
        is_system = (binary == "system")
        return (binary=binary, is_system_mpi=is_system)
    catch e
        # Fallback if MPIPreferences is not available
        return (binary="unknown", is_system_mpi=nothing)
    end
end

end # module