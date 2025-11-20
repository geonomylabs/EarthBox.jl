module MumpsChecks

import EarthBox.SecurityUtils: check_environment_variables
import EarthBox.SecurityUtils: check_file_permissions
import EarthBox.PathValidation: validate_directory_path
import EarthBox.PathValidation: validate_file_path
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState

function validate_environment_variables()::Nothing
    # Check required environment variables
    required_env_vars = [
        "JULIA_MUMPS_LIBRARY_PATH",
        "LD_LIBRARY_PATH"
    ]
    check_environment_variables(required_env_vars)
    return nothing
end

function validate_system_of_equations_directory(soe_dir_path::String)::Nothing
    validated_soe_dir = validate_directory_path(soe_dir_path)
    if !check_file_permissions(validated_soe_dir, "r")
        error("Insufficient permissions for system of equations directory: $validated_soe_dir")
    end
    return nothing
end

function validate_src_dir(src_dir_path::String)::String
    validated_src_dir = validate_directory_path(src_dir_path)
    if !check_file_permissions(validated_src_dir, "r")
        error("Insufficient permissions for source directory: $validated_src_dir")
    end
    return validated_src_dir
end

function validate_path_to_mumps_script(path_to_mumps_script::String)::String
    validate_file_path(path_to_mumps_script)
    if !check_file_permissions(path_to_mumps_script, "r")
        error("Insufficient permissions for MUMPS script: $path_to_mumps_script")
    end
    return path_to_mumps_script
end

function validate_earthbox_directory_from_src_dir(src_dir::String)::String
    path_to_earthbox = dirname(src_dir)
    validated_earthbox_dir = validate_directory_path(path_to_earthbox)
    if !check_file_permissions(validated_earthbox_dir, "r")
        error("Insufficient permissions for EarthBox directory: $validated_earthbox_dir")
    end
    return validated_earthbox_dir
end

function validate_solver_inputs(solver_config::SolverConfigState)::Nothing
    check_for_security_issues_with_solver_inputs(solver_config)
    check_nproc(solver_config.nproc)
    check_nthreads(solver_config.nthreads)
    check_verbose_output(solver_config.verbose_output)
    return nothing
end

function check_for_security_issues_with_solver_inputs(solver_config::SolverConfigState)
    check_integer_input(solver_config.nproc, "nproc")
    check_integer_input(solver_config.nthreads, "nthreads")
    check_analysis_method_input(solver_config.analysis_method)
    check_parallel_ordering_method_input(solver_config.parallel_ordering_method)
    check_integer_input(solver_config.verbose_output, "verbose_output")
    check_integer_or_float_input(solver_config.memory_relax_perc, "memory_relax_perc")
end

function check_nproc(nproc::Int)
    if nproc < 1
        error("Parameter nproc in solver.yml must be greater than 0.")
    end
end

function check_nthreads(nthreads::Int)
    if nthreads < 1
        error("Parameter nthreads in solver.yml must be greater than 0.")
    end
end

function check_verbose_output(verbose_output::Int)
    if verbose_output ∉ [0, 1]
        error("Parameter verbose_output in solver.yml must be 0 or 1.")
    end
end

function check_memory_relax_perc(memory_relax_perc::Int)
    if memory_relax_perc <= 0
        error("Parameter memory_relax_perc in solver.yml must be greater than 0.")
    end
end

function check_integer_or_float_input(param::Any, param_name::String)
    if !(isa(param, Int) || isa(param, Float64))
        error("Expected integer or float in solver.yml for $param_name but got $param")
    end
end

function check_integer_input(param::Any, param_name::String)
    if !isa(param, Int)
        error("Expected integer in solver.yml for $param_name but got $param")
    end
end

function check_analysis_method_input(param::String)
    if param ∉ ["PARALLEL", "SERIAL"]
        error("Expected PARALLEL or SERIAL in solver.yml but got $param")
    end
end

function check_parallel_ordering_method_input(param::String)
    if param ∉ ["PTSCOTCH", "ParMETIS"]
        error("Expected PTSCOTCH or ParMETIS in solver.yml but got $param")
    end
end

end