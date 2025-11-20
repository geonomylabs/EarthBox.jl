module RunMumpsSolverExternal

import Printf
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState
import EarthBox.PathValidation: validate_directory_path
import ..MumpsChecks

"""
    mpirun_mumps(output_dir::String, src_dir::String, solver_config::SolverConfigState):: String

Call MUMPS subprocess function.

This function checks for the existence of the system of equations directory
and runs the pymumps solver using the run_mumps_solver function.

# Arguments
- `output_dir::String`: Path to the output directory where the system of equations is stored
- `src_dir::String`: Path to the source directory containing the parallel_solver package
- `solver_config::SolverConfigState`: Solver configuration parameters

# Returns
- `String`: Standard output from running the MUMPS subprocess
"""
function mpirun_mumps(
    output_dir::String,
    src_dir::String,
    solver_config::SolverConfigState
)::String
    validate_directory_path(output_dir)
    validate_directory_path(src_dir)
    soe_dir_path = joinpath(output_dir, "system_of_equations")
    if isdir(soe_dir_path)
        output = run_mumps_solver(soe_dir_path, src_dir, solver_config)
        return output
    else
        error("System of equations directory not found.")
    end
end

"""
    run_mumps_solver(soe_dir_path::String, src_dir_path::String, solver_config::SolverConfigState) -> String

Run pymumps via mpirun using subprocess module.

# Arguments
- `soe_dir_path::String`: Path to directory containing npy files that define the system of equations
- `src_dir_path::String`: Path to source directory containing parallel_solver package
- `solver_config::SolverConfigState`: Solver configuration parameters

# Returns
- `String`: Standard output from running MUMPS subprocess
"""
function run_mumps_solver(
    soe_dir_path::String,
    src_dir_path::String,
    solver_config::SolverConfigState
)::String
    path_to_mumps_script = security_checks(soe_dir_path, src_dir_path, solver_config)
    cmd = make_command(solver_config, path_to_mumps_script, soe_dir_path, src_dir_path)
    output = execute_command(cmd, solver_config.pymumps_timeout)
    return output
end

function security_checks(
    soe_dir_path::String,
    src_dir_path::String,
    solver_config::SolverConfigState
)::String
    MumpsChecks.validate_system_of_equations_directory(soe_dir_path)
    path_to_mumps_script = define_path_to_external_mumps_solver_script(src_dir_path)
    MumpsChecks.validate_environment_variables()
    MumpsChecks.validate_solver_inputs(solver_config)
    return path_to_mumps_script
end

function define_path_to_external_mumps_solver_script(src_dir_path::String)::String
    validated_src_dir = MumpsChecks.validate_src_dir(src_dir_path)
    path_to_mumps_script = joinpath(
        validated_src_dir, 
        "parallel_solver", 
        "mumps_scripts", 
        "mumps_solver_external.jl"
        )
    path_to_mumps_script = MumpsChecks.validate_path_to_mumps_script(path_to_mumps_script)
    return path_to_mumps_script
end

function make_command(
    solver_config::SolverConfigState,
    path_to_mumps_script::String,
    soe_dir_path::String,
    src_dir_path::String
)::Cmd
    nproc = solver_config.nproc
    nthreads = solver_config.nthreads
    analysis_method = solver_config.analysis_method
    parallel_ordering_method = solver_config.parallel_ordering_method
    verbose_output = solver_config.verbose_output
    memory_relax_perc = solver_config.memory_relax_perc
    mpi_cmd = `mpirun -np $nproc`
    julia_cmd = `julia`
    script_path = `$path_to_mumps_script`
    args = [
        soe_dir_path,
        src_dir_path,
        string(nthreads),
        analysis_method,
        parallel_ordering_method,
        string(verbose_output),
        string(memory_relax_perc)
    ]
    cmd = `$mpi_cmd $julia_cmd $script_path $args`
    return cmd
end

function execute_command(cmd::Cmd, pymumps_timeout::Float64)::String
    try
        output_buffer = IOBuffer()
        process = run(pipeline(cmd, stdout=IOContext(
            output_buffer, :color=>true), stderr=stderr), wait=false)
        timed_out = false
        start_time = time()
        while process_running(process)
            delta_time = time() - start_time
            if delta_time > pymumps_timeout
                timed_out = true
                break
            end
            sleep(0.5)
        end
        if timed_out
            kill(process)
            Printf.@printf(
                "!!! WARNING !!! MUMPS.jl timed out after %d seconds.\n",
                pymumps_timeout
                )
            return "MUMPS.jl timed out."
        end
        seekstart(output_buffer)
        captured_output = String(take!(output_buffer))
        return captured_output
    catch e
        println(">> Error in command execution function. error: ", e)
        if isa(e, ProcessFailedException)
            Printf.@printf("!!! WARNING !!! MUMPS.jl failed: %s\n", e)
            return "MUMPS.jl failed."
        else
            rethrow(e)
        end
    end
end

end # module