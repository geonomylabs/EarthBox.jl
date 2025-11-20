module ConfigurationManager

include("configs/SolverConfig.jl")
include("configs/OutputConfig.jl")

import MPI
import .SolverConfig: SolverConfigState
import .OutputConfig: OutputConfigState

struct Configuration
    solver::SolverConfigState
    output::OutputConfigState
    function Configuration(
        run_paths::Union{Dict{String, String}, Nothing}, 
        use_mumps::Bool, 
        nprocs::Int,
        use_internal_mumps::Bool,
        mpi_comm::Union{MPI.Comm, Nothing},
        mpi_initialized::Bool,
        mpi_rank::Int,
        pass_large_arrays_via_mpi::Bool
    )
        if run_paths === nothing
            solver = SolverConfigState()
            output = OutputConfigState()
        else
            solver = SolverConfigState(
                use_mumps=use_mumps,
                nproc=nprocs,
                output_dir=run_paths["output_dir"],
                src_dir=run_paths["src_dir"],
                use_internal_mumps=use_internal_mumps,
                mpi_comm=mpi_comm,
                mpi_initialized=mpi_initialized,
                mpi_rank=mpi_rank,
                pass_large_arrays_via_mpi=pass_large_arrays_via_mpi
            )
            output = OutputConfigState()
        end
        return new(solver, output)
    end
end

end #module