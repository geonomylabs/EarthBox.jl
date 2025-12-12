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
    soe_dir_path = ARGS[1]
    pass_large_arrays_via_mpi = parse(Bool, ARGS[2])
    t1 = time()
    if pass_large_arrays_via_mpi
        solve_system_mpi_comm(soe_dir_path)
    else
        solve_system_io_comm(soe_dir_path)
    end
    t2 = time()
    println("Time to solve system:", t2 - t1)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    MumpsSolverLoop.main()
end