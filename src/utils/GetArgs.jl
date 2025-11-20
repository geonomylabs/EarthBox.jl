module GetArgs

using TOML

export get_model_output_path_from_cl_args

"""
    get_model_output_path_from_cl_args(model_output_path_from_script_path::String)::String

Get the model output path from the command line argument `model_output_path=...`.
If the model output path is not provided, then the model output path is set to the 
path from the script executing the plotter.

# Arguments
- `model_output_path_from_script_path::String`:
    - Path to the model output directory sent from the script executing the plotter.

# Returns
- `model_output_path::String`:
    - Path to the model output directory

"""
function get_model_output_path_from_cl_args(model_output_path_from_script_path::String)::String
    model_output_path = get_model_output_path_from_args()
    if isnothing(model_output_path)
        model_output_path = model_output_path_from_script_path
    end
    return model_output_path
end

"""
    get_case_name_from_cl_args()::String

Get the model case name from the command line argument `case_name=...`.
If the model case name is not provided, then the model case name is set to "case0".

# Returns
- `model_case_name::String`:
    - Name of the model case

"""
function get_case_name_from_cl_args()::String
    model_case_name = get_model_case_name_from_args()
    if isnothing(model_case_name)
        model_case_name = "case0"
    end
    return model_case_name
end

"""
    get_istart()::Int64

Get the starting time step from the command line argument `istart=...`.
If the starting time step is not provided, then the starting time step is set to 1.

# Returns
- `istart::Int64`:
    - Starting time step

"""
function get_istart()::Int64
    istart = get_istart_from_args()
    if isnothing(istart)
        istart = 1
    end
    return istart
end

"""
    get_iend()::Union{Int64, Nothing}

Get the ending time step from the command line argument `iend=...`.
If the ending time step is not provided, then the ending time step is set to nothing.

# Returns
- `iend::Union{Int64, Nothing}`:
    - Ending time step

"""
function get_iend()::Union{Int64, Nothing}
    iend = get_iend_from_args()
    if isnothing(iend)
        iend = nothing
    end
    return iend
end

"""
    get_plot_option_name()::String

Get the standardized plot option name from the command line argument.
If the plot option name is not provided, then an error is thrown.

The plot option name is used to lookup a user defined function provided by the 
user for plotting EarthBox output.

# Valid Plot Option Names
- `marker_plots`: Plot markers using user defined function.
- `scalar_plots`: Plot scalars using user defined function.
- `velocity_plots`: Plot velocity using user defined function.
- `stokes_convergence_plots`: Plot Stokes convergence using user defined function.
- `yield_strength_plots`: Plot yield strength using user defined function.

# Returns
- `plot_option_name::String`:
    - Name of the standardized plot option

"""
function get_plot_option_name()::String
    plot_option_name = nothing
    if "marker_plots" in ARGS
        plot_option_name = "marker_plots"
    end
    if "scalar_plots" in ARGS
        plot_option_name = "scalar_plots"
    end
    if "stokes_convergence_plots" in ARGS
        plot_option_name = "stokes_convergence_plots"
    end
    if "velocity_plots" in ARGS
        plot_option_name = "velocity_plots"
    end
    if "yield_strength_plots" in ARGS
        plot_option_name = "yield_strength_plots"
    end
    if isnothing(plot_option_name)
        throw(ArgumentError("No plot option name provided"))
    end
    if !(plot_option_name in ["marker_plots", "scalar_plots", "stokes_convergence_plots", "velocity_plots", "yield_strength_plots"])
        throw(ArgumentError("$(plot_option_name) is not a valid option."))
    end
    return plot_option_name
end

"""
    manage_earthbox_project_path(
        eb_project_path_from_script::Union{String, Nothing} = nothing
    )::Union{String, Nothing}

Get the EarthBox project path from command line arguments or default path
if provided by the script executing the action.

# Returns
- `eb_path::Union{String, Nothing}`:
    - Path to the EarthBox project

# Arguments
- `eb_project_path_from_script::Union{String, Nothing}`:
    - Path to the EarthBox project from the script executing the action.

"""
function manage_earthbox_project_path(
    eb_project_path_from_script::Union{String, Nothing} = nothing
)::Union{String, Nothing}
    eb_path = get_earthbox_project_path_from_args()
    if !isnothing(eb_project_path_from_script) && !isnothing(eb_path)
        return eb_project_path_from_script
    end
    return eb_path
end

"""
    get_earthbox_project_path_from_args()::Union{String, Nothing}

Get the EarthBox project path from the command line argument `eb_path=...`.
If the EarthBox project path is not provided, then the EarthBox project path is 
set to nothing.

# Returns
- `eb_path::Union{String, Nothing}`:
    - Path to the EarthBox project

"""
function get_earthbox_project_path_from_args()::Union{String, Nothing}
    for arg in ARGS
        if startswith(arg, "eb_path=") || startswith(arg, "eb_path = ") || startswith(arg, "eb_path = ") || startswith(arg, "eb_path= ")
            # Remove "eb_path=" or "eb_path = " prefix and any leading/trailing whitespace
            return strip(replace(arg, r"^eb_path\s*=\s*" => ""))
        end
    end
    return nothing
end

function get_model_case_name()::String
    if length(ARGS) > 0 && ARGS[1] isa String
        return ARGS[1]
    else
        return "case0"
    end
end

function get_model_case_name_for_plotting()::String
    if length(ARGS) < 2
        return "case0"
    else
        return ARGS[2]
    end
end

function get_runit_actions_from_cl_args()
    run_model = false
    plot_markers = false
    plot_scalars = false
    if "run_model" in ARGS
        run_model = true
    end
    if "plot_markers" in ARGS
        plot_markers = true
    end
    if "plot_scalars" in ARGS
        plot_scalars = true
    end
    return run_model, plot_markers, plot_scalars
end

function get_root_path_from_args()::Union{String, Nothing}
    for arg in ARGS
        if startswith(arg, "root_path=") || startswith(arg, "root_path = ") || startswith(arg, "root_path = ") || startswith(arg, "root_path= ")
            # Remove "root_path=" or "root_path = " prefix and any leading/trailing whitespace
            return strip(replace(arg, r"^root_path\s*=\s*" => ""))
        end
    end
    return nothing
end

function get_storage_path_from_args()::Union{String, Nothing}
    for arg in ARGS
        if startswith(arg, "storage_path=") || startswith(arg, "storage_path = ") || startswith(arg, "storage_path = ") || startswith(arg, "storage_path= ")
            # Remove "storage_path=" or "storage_path = " prefix and any leading/trailing whitespace
            return strip(replace(arg, r"^storage_path\s*=\s*" => ""))
        end
    end
    return nothing
end

function get_model_output_path_from_args()::Union{String, Nothing}
    for arg in ARGS
        if startswith(arg, "model_output_path=") || startswith(arg, "model_output_path = ") || startswith(arg, "model_output_path = ") || startswith(arg, "model_output_path= ")
            # Remove "model_output_path=" or "model_output_path = " prefix and any leading/trailing whitespace
            return strip(replace(arg, r"^model_output_path\s*=\s*" => ""))
        end
    end
    return nothing
end

function get_model_output_path_from_args_strict()::String
    for arg in ARGS
        if startswith(arg, "model_output_path=") || startswith(arg, "model_output_path = ") || startswith(arg, "model_output_path = ") || startswith(arg, "model_output_path= ")
            # Remove "model_output_path=" or "model_output_path = " prefix and any leading/trailing whitespace
            return strip(replace(arg, r"^model_output_path\s*=\s*" => ""))
        end
    end
    throw(ArgumentError("Model output path not found in command line arguments"))
end

function get_model_case_name_from_args()::Union{String, Nothing}
    for arg in ARGS
        if startswith(arg, "case_name=") || startswith(arg, "case_name = ") || startswith(arg, "case_name = ") || startswith(arg, "case_name= ")
            # Remove "case_name=" or "case_name = " prefix and any leading/trailing whitespace
            return strip(replace(arg, r"^case_name\s*=\s*" => ""))
        end
    end
    return nothing
end

function get_istart_from_args()::Union{Int64, Nothing}
    for arg in ARGS
        if startswith(arg, "istart=") || startswith(arg, "istart = ") || startswith(arg, "istart = ") || startswith(arg, "istart= ")
            # Remove "istart=" or "istart = " prefix and any leading/trailing whitespace
            return parse(Int64, strip(replace(arg, r"^istart\s*=\s*" => "")))
        end
    end
    return nothing
end

function get_iend_from_args()::Union{Int64, Nothing}
    for arg in ARGS
        if startswith(arg, "iend=") || startswith(arg, "iend = ") || startswith(arg, "iend = ") || startswith(arg, "iend= ")
            # Remove "iend=" or "iend = " prefix and any leading/trailing whitespace
            return parse(Int64, strip(replace(arg, r"^iend\s*=\s*" => "")))
        end
    end
    return nothing
end

function get_model_input_file_from_args()::String
    for arg in ARGS
        if startswith(arg, "model_input_file=") || startswith(arg, "model_input_file = ") || startswith(arg, "model_input_file = ") || startswith(arg, "model_input_file= ")
            # Remove "model_input_file=" or "model_input_file = " prefix and any leading/trailing whitespace
            path = strip(replace(arg, r"^model_input_file\s*=\s*" => ""))
            if isfile(path)
                # Check if the file has ".yml" extension (lowercase only)
                if endswith(path, ".yml")
                    return path
                end
                throw(ArgumentError("Model input file must have .yml extension: $(path)"))
            end
            throw(ArgumentError("Model input file not found: $(path)"))
        end
    end
    throw(ArgumentError("Model input file not found in command line arguments"))
end

function get_materials_input_file_from_args()::String
    for arg in ARGS
        if startswith(arg, "materials_input_file=") || startswith(arg, "materials_input_file = ") || startswith(arg, "materials_input_file = ") || startswith(arg, "materials_input_file= ")
            # Remove "materials_input_file=" or "materials_input_file = " prefix and any leading/trailing whitespace
            path = strip(replace(arg, r"^materials_input_file\s*=\s*" => ""))
            if isfile(path)
                # Check if the file has ".yml" extension (lowercase only)
                if endswith(path, ".yml")
                    return path
                end
                throw(ArgumentError("Materials input file must have .yml extension: $(path)"))
            end
            throw(ArgumentError("Materials input file not found: $(path)"))
        end
    end
    throw(ArgumentError("Materials input file not found in command line arguments"))
end

function get_materials_library_file_from_args()::String
    for arg in ARGS
        if startswith(arg, "materials_library_file=") || startswith(arg, "materials_library_file = ") || startswith(arg, "materials_library_file = ") || startswith(arg, "materials_library_file= ")
            # Remove "materials_library_file=" or "materials_library_file = " prefix and any leading/trailing whitespace
            path = strip(replace(arg, r"^materials_library_file\s*=\s*" => ""))
            if isfile(path)
                # Check if the file has ".yml" extension (lowercase only)
                if endswith(path, ".yml")
                    return path
                end
                throw(ArgumentError("Materials library file must have .yml extension: $(path)"))
            end
            throw(ArgumentError("Materials library file not found: $(path)"))
        end
    end
    throw(ArgumentError("Materials library file not found in command line arguments"))
end

function get_make_backup_from_args()::Bool
    for arg in ARGS
        if startswith(arg, "make_backup=") || startswith(arg, "make_backup = ") || startswith(arg, "make_backup = ") || startswith(arg, "make_backup= ")
            # Remove "make_backup=" or "make_backup = " prefix and any leading/trailing whitespace
            return parse(Bool, strip(replace(arg, r"^make_backup\s*=\s*" => "")))
        end
    end
    return false
end

function get_restart_from_backup_from_args()::Bool
    for arg in ARGS
        if startswith(arg, "restart_from_backup=") || startswith(arg, "restart_from_backup = ") || startswith(arg, "restart_from_backup = ") || startswith(arg, "restart_from_backup= ")
            # Remove "restart_from_backup=" or "restart_from_backup = " prefix and any leading/trailing whitespace
            return parse(Bool, strip(replace(arg, r"^restart_from_backup\s*=\s*" => "")))
        end
    end
    return false
end

function get_use_mumps_from_args()::Bool
    for arg in ARGS
        if startswith(arg, "use_mumps=") || startswith(arg, "use_mumps = ") || startswith(arg, "use_mumps = ") || startswith(arg, "use_mumps= ")
            # Remove "use_mumps=" or "use_mumps = " prefix and any leading/trailing whitespace
            return parse(Bool, strip(replace(arg, r"^use_mumps\s*=\s*" => "")))
        end
    end
    return false
end

function get_use_mumps_internal_from_args()::Bool
    for arg in ARGS
        if startswith(arg, "use_mumps_internal=") || startswith(arg, "use_mumps_internal = ") || startswith(arg, "use_mumps_internal = ") || startswith(arg, "use_mumps_internal= ")
            # Remove "use_mumps_internal=" or "use_mumps_internal = " prefix and any leading/trailing whitespace
            return parse(Bool, strip(replace(arg, r"^use_mumps_internal\s*=\s*" => "")))
        end
    end
    return false
end

function get_nprocs_from_args()::Int
    for arg in ARGS
        if startswith(arg, "nprocs=") || startswith(arg, "nprocs = ") || startswith(arg, "nprocs = ") || startswith(arg, "nprocs= ")
            # Remove "nprocs=" or "nprocs = " prefix and any leading/trailing whitespace
            return parse(Int, strip(replace(arg, r"^nprocs\s*=\s*" => "")))
        end
    end
    return 1
end

function get_marker_output_dict_from_args()::Union{Dict{String, Bool}, Nothing}
    raw_dict = get_arg_value("marker_output")
    if isnothing(raw_dict)
        return nothing
    end
    marker_output_dict = parse_inline_toml_dict(raw_dict)
    return marker_output_dict
end

function get_arg_value(name::String)::Union{String, Nothing}
    raw_dict = nothing
    for (i, arg) in pairs(ARGS)
        if startswith(arg, name * "=")
            raw_dict = split(arg, "=", limit=2)[2]
        elseif arg == name
            if i < length(ARGS)
                next_arg = ARGS[i+1]
                if startswith(next_arg, "=")
                    if i+1 < length(ARGS)
                        next_next_arg = ARGS[i+2]
                        if startswith(next_next_arg, "{")
                            raw_dict = next_next_arg
                        else
                            raw_dict = nothing
                        end
                    else
                        raw_dict = nothing
                    end
                else
                    raw_dict = next_arg
                end
            else
                raw_dict = nothing
            end
        end
    end
    return raw_dict
end

function parse_inline_toml_dict(raw_dict::String)::Union{Dict{String, Bool}, Nothing}
    doc = Dict{String, Any}()
    try
        doc = TOML.parse("marker_output = " * raw_dict)
    catch e
        println("!!! WARNING !!! Parsing inline TOML dictionary failed. Type: $(typeof(e)) with message: $(e)")
        println("Make sure there are no spaces between the marker_output argument name or value and = sign.")
        return nothing
    end
    dict = doc["marker_output"]
    dict_final = Dict{String, Bool}()
    for (key, value) in dict
        if value isa Bool
            dict_final[key] = Bool(value)
        elseif value isa String && value in ["true", "false"]
            dict_final[key] = parse(Bool, value)
        end
    end
    return dict_final
end

end # module GetArgs