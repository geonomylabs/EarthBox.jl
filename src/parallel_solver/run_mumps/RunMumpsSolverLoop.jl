module RunMumpsSolverLoop

import Printf
import MPI
import EarthBox.MPIManager: get_mpi_implementation
import EarthBox.PrintFuncs: print_info, print_warning
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState
import EarthBox.ConfigurationManager.SolverConfig: InternalMumpsSolver
import LinearAlgebra: norm
import ..MumpsChecks
import ..NamesManager

function initialize_persistent_solver(
    solver_config::SolverConfigState,
    comm::Union{MPI.Comm, Nothing}
)::Bool
    solver = solver_config.internal_mumps_solver
    use_mumps = solver.use_mumps

    mpi_implementation = get_mpi_implementation()
    is_system_mpi = mpi_implementation.is_system_mpi
    print_info("MPI implementation is system?: $(is_system_mpi)", level=1)
    if !is_system_mpi && use_mumps
        # Need to add a check for the mpiexecjl command
    end

    if use_mumps && comm !== nothing
        # Need to remove pre-existing solver flag files that send signals to the 
        # solver loop to start solving or terminate
        clean_up_soe_dir(solver_config)
        if solver.is_running
            return true  # Already initialized
        end
        try
            if is_system_mpi
                solver.inter = run_earthbox_mumps_solver_using_system_mpi(solver_config, comm)
            else
                run_earthbox_mumps_solver_using_mpiexecjl(solver_config)
            end
            solver.is_running = true
            return true
        catch e
            print_warning("Failed to initialize persistent MUMPS solver: $e")
            solver.is_running = false
            return false
        end
    end
    return false
end

function restart_persistent_solver(
    solver_config::SolverConfigState,
    comm::Union{MPI.Comm, Nothing}
)::Bool
    solver = solver_config.internal_mumps_solver
    solver.is_running = false
    check = initialize_persistent_solver(solver_config, comm)
    return check
end

function run_earthbox_mumps_solver_using_system_mpi(
    solver_config::SolverConfigState,
    comm::MPI.Comm
)::MPI.Comm
    path_to_mumps_script = security_checks(solver_config)
    args = get_args(solver_config, path_to_mumps_script)
    inter = spawn_external_mumps_solver_using_system_mpi(args, solver_config, comm)
    return inter
end

function run_earthbox_mumps_solver_using_mpiexecjl(
    solver_config::SolverConfigState
)::String
    path_to_mumps_script = security_checks(solver_config)
    mpi_project_path = nothing
    # Only set mpi_project_path if MPI.jl is not available in Julia's shared environment
    # Check if MPI is in the default/shared environment (e.g., ~/.julia/environments/v1.9/)
    is_base_mpi = is_mpi_in_shared_environment()
    print_info("Is MPI.jl found in the base Julia environment?: $(is_base_mpi)", level=1)
    if !is_base_mpi
        # MPI.jl is not in shared environment, so we need to specify the project path
        mpi_project_path = dirname(solver_config.src_dir)
        print_info("MPI.jl not found in shared environment, using project path: $(mpi_project_path)", level=2)
    else
        print_info("MPI.jl found in shared environment, no project path needed", level=2)
    end
    output = spawn_external_mumps_solver_using_mpiexecjl(
        solver_config, path_to_mumps_script, mpi_project_path)
    return output
end

"""

    is_mpi_in_shared_environment()::Bool

Check if MPI.jl is available in Julia's shared/default environment.
Returns true if MPI is directly available in the shared environment, false otherwise.
Note: We check Project.toml because Project.toml lists packages.

"""
function is_mpi_in_shared_environment()::Bool
    try
        # Get Julia version (e.g., "1.9")
        version_string = "v$(VERSION.major).$(VERSION.minor)"
        print_info("Julia version: $(version_string)", level=1)
        # Construct path to shared environment Project.toml (not Manifest.toml)
        shared_env_path = joinpath(DEPOT_PATH[1], "environments", version_string, "Project.toml")
        print_info("Shared environment project path: $(shared_env_path)", level=1)
        # Check if Project.toml file exists
        if !isfile(shared_env_path)
            print_info("Shared environment Project.toml not found", level=1)
            return false
        end
        # Read the Project.toml file and check if MPI is listed as a direct dependency
        project_content = read(shared_env_path, String)
        # Look for MPI in the [deps] section
        # Format is: MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
        check = occursin(r"(?:^|\n)MPI\s*=", project_content)
        print_info("MPI found in shared environment Project.toml?: $(check)", level=1)
        return check
    catch e
        # If any error occurs, assume MPI is not in shared environment
        print_warning("Error checking shared environment for MPI: $e", level=2)
        return false
    end
end

function spawn_external_mumps_solver_using_system_mpi(
    args::Vector{String}, 
    solver_config::SolverConfigState,
    comm::MPI.Comm,
)::Union{MPI.Comm, String}
    N = solver_config.nproc
    try
        print_info("Spawning external MUMPS solver loop ...", level=1)
        # Spawn the MPI processes
        inter = MPI.Comm_spawn("julia", args, N, comm)
        print_info("External MUMPS solver loop spawned", level=1)
        return inter
    catch e
        if isa(e, ProcessFailedException)
            print_warning("MUMPS system_solver loop failed to spawn: $e")
            return "MUMPS failed to spawn: $(e)"
        else
            rethrow(e)
        end
    end
end

function spawn_external_mumps_solver_using_mpiexecjl(
    solver_config::SolverConfigState,
    path_to_mumps_scripts::String,
    mpi_project_path::Union{String, Nothing} = nothing
)::String
    N = solver_config.nproc
    try
        print_info("Spawning external MUMPS solver loop using mpiexecjl...", level=1)
        cmd = make_command(solver_config, path_to_mumps_scripts, mpi_project_path)
        println("Command: ", cmd)
        output = execute_command(cmd)
        print_info("External MUMPS solver loop spawned using mpiexecjl", level=1)
        return output
    catch e
        if isa(e, ProcessFailedException)
            print_warning("MUMPS mpiexecjl loop failed to spawn: $e")
            return "MUMPS failed to spawn: $(e)"
        else
            rethrow(e)
        end
    end
end

function make_command(
    solver_config::SolverConfigState,
    path_to_mumps_script::String,
    mpi_project_path::Union{String, Nothing}
)::Cmd
    nproc = solver_config.nproc

    analysis_method = solver_config.analysis_method
    parallel_ordering_method = solver_config.parallel_ordering_method
    verbose_output = solver_config.verbose_output
    memory_relax_perc = solver_config.memory_relax_perc
    pass_large_arrays_via_mpi = false
    # path_to_earthbox = dirname(solver_config.src_dir)
    
    no_startup = "--startup-file=no"
    if mpi_project_path !== nothing
        mpi_project = "--project=$mpi_project_path"
    else
        mpi_project = nothing
    end
    
    script_path = string(path_to_mumps_script)
    soe_dir_path = get_soe_dir_path(solver_config)

    args = [
        string(analysis_method),
        string(parallel_ordering_method),
        string(verbose_output),
        string(memory_relax_perc),
        soe_dir_path,
        string(pass_large_arrays_via_mpi),
    ]

    if mpi_project !== nothing
        cmd = `mpiexecjl $mpi_project -np $nproc julia $no_startup $script_path $args`
    else
        cmd = `mpiexecjl -np $nproc julia $no_startup $script_path $args`
    end
    
    # We now require that EarthBox is installed in a standard Julia location or 
    # the JULIA_PROJECT environment variable is set equal to the EarthBox project path.
    # In this case the `--project=<eb_path>` argument is not included in the command.
    # eb_project = "--project=$path_to_earthbox"
    # cmd = `mpiexecjl $mpi_project -np $nproc julia $eb_project $no_startup $script_path $args`
    
    return cmd
end

function execute_command(cmd::Cmd)::String
    try
        output_buffer = IOBuffer()
        process = run(pipeline(cmd, stdout=IOContext(output_buffer, :color=>true), stderr=stderr), wait=false)
        return "Spawned"
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

function security_checks(
    solver_config::SolverConfigState
)::String
    path_to_mumps_script = define_path_to_internal_mumps_solver_script(solver_config)
    MumpsChecks.validate_environment_variables()
    MumpsChecks.validate_solver_inputs(solver_config)
    return path_to_mumps_script
end

function define_path_to_internal_mumps_solver_script(solver_config::SolverConfigState)::String
    src_dir = solver_config.src_dir
    validated_earthbox_dir = MumpsChecks.validate_earthbox_directory_from_src_dir(src_dir)
    path_to_mumps_script = joinpath(
        validated_earthbox_dir, 
        "src", 
        "parallel_solver", 
        "mumps_scripts", 
        "mumps_solver_loop.jl"
    )
    path_to_mumps_script = MumpsChecks.validate_path_to_mumps_script(path_to_mumps_script)
    return path_to_mumps_script
end

function get_args(
    solver_config::SolverConfigState,
    path_to_mumps_script::String
)::Vector{String}
    analysis_method = solver_config.analysis_method
    parallel_ordering_method = solver_config.parallel_ordering_method
    verbose_output = solver_config.verbose_output
    memory_relax_perc = solver_config.memory_relax_perc
    pass_large_arrays_via_mpi = solver_config.pass_large_arrays_via_mpi

    path_to_earthbox = dirname(solver_config.src_dir)
    no_startup = "--startup-file=no"
    no_banner = "--banner=no"
    project = "--project=$path_to_earthbox"
    script_path = path_to_mumps_script
    soe_dir_path = get_soe_dir_path(solver_config)
    args = [
        project,
        no_startup,
        no_banner,
        script_path,
        analysis_method,
        parallel_ordering_method,
        string(verbose_output),
        string(memory_relax_perc),
        soe_dir_path,
        string(pass_large_arrays_via_mpi)
    ]
    return args
end

function mpirun_mumps_mpi_comm(
    solver_config::SolverConfigState,
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    rhs::Vector{Float64},
    matrix_info::Vector{Int64}
)::Tuple{Vector{Float64}, String}
    solver = solver_config.internal_mumps_solver
    if !solver.is_running
        error("Persistent solver is not running. Call initialize_persistent_solver() first.")
    end
    result = execute_mumps_with_timeout_mpi_comm(
        solver.inter, solver_config, Li, Lj, Lv, rhs, matrix_info)
    mumps_solver_status = failure_check(result)
    if mumps_solver_status == "Failure"
        S = zeros(Float64, 1)
    else
        S = result
    end
    return S, mumps_solver_status
end

function mpirun_mumps_io_comm(
    solver_config::SolverConfigState,
    termination_flag::Int64
)::String
    solver = solver_config.internal_mumps_solver
    if !solver.is_running
        error("Persistent solver is not running. Call initialize_persistent_solver() first.")
    end
    solution_flag = execute_mumps_with_timeout_io_comm(solver_config, termination_flag)
    if solution_flag == 0
        return "Failure"
    end
    return "Success"
end

# Shutdown the persistent solver
function shutdown_persistent_solver(solver_config::SolverConfigState)
    pass_large_arrays_via_mpi = solver_config.pass_large_arrays_via_mpi
    solver = solver_config.internal_mumps_solver
    if solver.is_running && solver.inter !== nothing
        try
            if MPI.Initialized()
                if pass_large_arrays_via_mpi
                    # Send termination signal (negative N)
                    termination_signal = [-1, 0]
                    MPI.Send(termination_signal, 0, 0, solver.inter)
                else
                    # termination communicated via file IO
                    soe_dir_path = get_soe_dir_path(solver_config)
                    check_system_of_equations_dir(soe_dir_path)
                    termination_info = -1
                    create_termination_file(soe_dir_path, termination_info)
                end
                solver.is_running = false
                print_info("Persistent MUMPS solver shutdown signal sent.")
            else
                print_warning("MPI already finalized, cannot send shutdown signal")
            end
            solver.is_running = false
        catch e
            print_warning("Error sending shutdown signal: $e", level=2)
        end
    end
end

# Check if solver is healthy
function is_solver_healthy(solver::InternalMumpsSolver)::Bool
    return solver.is_running && solver.inter !== nothing
end

function execute_mumps_with_timeout_io_comm(
    solver_config::SolverConfigState, 
    termination_flag::Int64
)::Union{Int64, String}
    print_info("Executing MUMPS solver with file IO communication...", level=2)
    pymumps_timeout = solver_config.pymumps_timeout
    loop_check_time = 1.0
    try
        soe_dir_path = get_soe_dir_path(solver_config)
        check_system_of_equations_dir(soe_dir_path)
        # Create a termination file to signal to the child process to terminate
        # solver loop
        create_termination_file(soe_dir_path, termination_flag)
        # Create a solver ready file to signal to the child process that it is time
        # to start solving the system of equations
        create_solver_ready_file(soe_dir_path)

        # Initialize solution flag to 0 (system not solved)
        solution_flag = 0

        timed_out = false
        start_time = time()
        last_print_time = start_time
        while !timed_out
            delta_print_time = time() - last_print_time
            delta_time = time() - start_time
            if delta_print_time >= loop_check_time
                print_info("Execute mumps loop active for $delta_time seconds...", level=2)
                last_print_time = time()
            end
            if delta_time > pymumps_timeout
                timed_out = true
                break
            end
            # Check for a solution flag file was produced by the external solver loop indicating 
            # that the solver is done solving the system of equations
            solution_flag_file_path = get_solution_flag_file_path(soe_dir_path)
            if isfile(solution_flag_file_path)
                # Read the solution flag from the binary solution flag file
                # 0 = solution was not produced possibly due to an error in the solver loop
                # 1 = solution was successfully produced and exported
                solution_flag = get_solution_flag_from_file(solution_flag_file_path)
                # Remove the solution flag file since it has been read
                rm(solution_flag_file_path)
                # Return the solution flag (1 if solution was found, 0 if not)
                return solution_flag
            end
            sleep(0.01)
        end
        if timed_out
            print_warning("MUMPS system_solver timed out after $pymumps_timeout seconds.", level=2)
            return "MUMPS timed out"
        end
    catch e
        if isa(e, ProcessFailedException)
            print_warning("MUMPS system_solver failed: $e", level=2)
            return "MUMPS failed: $(e)"
        else
            rethrow(e)
        end
    end
end

function check_system_of_equations_dir(soe_dir_path::String)::Nothing
    if !isdir(soe_dir_path)
        mkpath(soe_dir_path)
        print_info("Created directory: $soe_dir_path", level=2)
    end
    return nothing
end

function get_soe_dir_path(solver_config::SolverConfigState)::String
    output_dir = solver_config.output_dir
    names = NamesManager.FileAndDirNames()
    soe_dir_path = joinpath(output_dir, names.soe_dir_name)
    return soe_dir_path
end

function get_solution_flag_file_path(soe_dir_path::String)::String
    names = NamesManager.FileAndDirNames()
    solution_file_name = names.solution_flag_file_name
    solution_file_path = joinpath(soe_dir_path, solution_file_name)
    return solution_file_path
end

function create_termination_file(
    soe_dir_path::String,
    termination_flag::Int64
)::Nothing
    names = NamesManager.FileAndDirNames()
    termination_file = joinpath(soe_dir_path, names.termination_info_file_name)
    open(termination_file, "w") do f
        write(f, termination_flag)
    end
    return nothing
end

function create_solver_ready_file(soe_dir_path::String)::Nothing
    names = NamesManager.FileAndDirNames()
    flag_file = joinpath(soe_dir_path, names.ready_to_solve_file_name)
    open(flag_file, "w") do f
        write(f, "ready")
    end
    return nothing
end

function clean_up_soe_dir(solver_config::SolverConfigState)::Nothing
    soe_dir_path = get_soe_dir_path(solver_config)
    remove_termination_file(soe_dir_path)
    remove_ready_to_solve_file(soe_dir_path)
    remove_solution_flag_file(soe_dir_path)
    return nothing
end

function remove_termination_file(soe_dir_path::String)::Nothing
    names = NamesManager.FileAndDirNames()
    file_name = names.termination_info_file_name
    file_path = joinpath(soe_dir_path, file_name)
    if isfile(file_path)
        rm(file_path)
    end
    return nothing
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

function remove_solution_flag_file(soe_dir_path::String)::Nothing
    names = NamesManager.FileAndDirNames()
    file_name = names.solution_flag_file_name
    file_path = joinpath(soe_dir_path, file_name)
    if isfile(file_path)
        rm(file_path)
    end
    return nothing
end

function get_solution_flag_from_file(solution_flag_file_path::String)::Int64
    if solution_flag_file_has_content(solution_flag_file_path)
        try
            open(solution_flag_file_path, "r") do f
                solution_flag = read(f, Int64)
                return solution_flag
            end
        catch e
            print_warning("Error reading solution flag from file: $e", level=2)
            return 0
        end
    else
        print_warning("Solution flag file is not present or empty: $solution_flag_file_path", level=2)
        return 0
    end
end

function solution_flag_file_has_content(solution_flag_file_path)::Bool
    if isfile(solution_flag_file_path)
        if filesize(solution_flag_file_path) >= sizeof(Int64)
            return true
        else
            return false
        end
    else
        return false
    end
    return false
end

function execute_mumps_with_timeout_mpi_comm(
    inter::MPI.Comm,
    solver_config::SolverConfigState, 
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    rhs::Vector{Float64},
    matrix_info::Vector{Int64}
)
    print_info("Executing MUMPS solver with timeout...", level=2)
    pymumps_timeout = solver_config.pymumps_timeout
    loop_check_time = 1.0
    try
        #MPI.Send(matrix_info, 0, 0, inter)
        req = MPI.Isend(matrix_info, 0, 0, inter)
        MPI.Wait(req)
        MPI.Send(Li, 0, 1, inter)
        MPI.Send(Lj, 0, 2, inter)
        MPI.Send(Lv, 0, 3, inter)
        MPI.Send(rhs, 0, 4, inter)

        N = matrix_info[1]
        x = Vector{Float64}(undef, N)

        timed_out = false
        start_time = time()
        last_print_time = start_time
        while !timed_out
            delta_print_time = time() - last_print_time
            delta_time = time() - start_time
            if delta_print_time >= loop_check_time
                print_info("Execute mumps loop active for $delta_time seconds...", level=2)
                last_print_time = time()
            end
            if delta_time > pymumps_timeout
                timed_out = true
                break
            end
            try
                has_solution, status = MPI.Iprobe(0, 5, inter)
                if has_solution
                    MPI.Recv!(x, 0, 5, inter)
                    if !solution_failed(x)
                        print_info("Solution min: $(minimum(x)), max: $(maximum(x)), norm: $(norm(x))", level=2)
                    end
                    return x
                end
            catch e
                print_warning("MPI communication of solution from solve_system failed: $e", level=2)
                return "MUMPS failed: $(e)"
            end
            sleep(0.01)  # Check every 100ms
        end
        if timed_out
            print_warning("MUMPS system_solver timed out after $pymumps_timeout seconds.", level=2)
            return "MUMPS timed out"
        end
    catch e
        if isa(e, ProcessFailedException)
            print_warning("MUMPS system_solver failed: $e", level=2)
            return "MUMPS failed: $(e)"
        else
            rethrow(e)
        end
    end
end

function failure_check(result)::String
    failure_status = "Success"
    if is_failure(result)
        print_warning("MUMPS solver failed or timed out: $result", level=2)
        failure_status = "Failure"
    end
    return failure_status
end

function is_failure(result)
    if isa(result, String) && (occursin("failed", result) || occursin("timed out", result) || occursin("MUMPS Error", result))
        return true
    end
    if isa(result, Vector{Float64}) && solution_failed(result)
        return true
    end
    return false
end

function is_success(result)
    if isa(result, Vector{Float64})
        return true
    end
    return false
end

function solution_failed(S::Vector{Float64})
    if isnan(S[1])
        return true
    end
    return false
end


end # module