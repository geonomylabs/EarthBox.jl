module MumpsSolverLoop

try
    old_stderr = stderr
    redirect_stderr(devnull)
    import EarthBox.ParallelSolver.MumpsLoop: get_integer_option_values
    import EarthBox.ParallelSolver.MumpsLoop: solve_system_io_comm
    import EarthBox.ParallelSolver.MumpsLoop: solve_system_mpi_comm
    redirect_stderr(old_stderr)
catch e
    println("!!! ERROR !!! Error during module import in mumps_solver_loop.jl. Comment out message suppression and re-run.")
    redirect_stderr(old_stderr)
    rethrow(e)
end

function main()
    analysis_method = ARGS[1]
    parallel_ordering_method = ARGS[2]
    verbose_output_itype = parse(Int, ARGS[3])
    memory_relax_perc = parse(Int, ARGS[4])
    soe_dir_path = ARGS[5]
    pass_large_arrays_via_mpi = parse(Bool, ARGS[6])
    (
        analysis_method_itype, 
        parallel_ordering_method_itype
    ) = get_integer_option_values(
        analysis_method,
        parallel_ordering_method
    )
    t1 = time()
    if pass_large_arrays_via_mpi
        solve_system_mpi_comm(
            analysis_method_itype,
            parallel_ordering_method_itype,
            verbose_output_itype,
            memory_relax_perc,
            soe_dir_path
        )
    else
        solve_system_io_comm(
            analysis_method_itype,
            parallel_ordering_method_itype,
            verbose_output_itype,
            memory_relax_perc,
            soe_dir_path
        )
    end
    t2 = time()
    println("Time to solve system:", t2 - t1)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    MumpsSolverLoop.main()
end