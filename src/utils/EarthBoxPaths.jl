module EarthBoxPaths

using Printf
import EarthBox.PrintFuncs: print_info

struct EarthBoxPathNames
    model_input_file::String
    materials_input_file::String
    materials_library_file::String
    output_dir::String

    EarthBoxPathNames() = new(
        "model_input_file",
        "materials_input_file",
        "materials_library_file",
        "output_dir"
    )
end

mutable struct EarthBoxPathsState
    key_names::EarthBoxPathNames
    required_paths::Vector{String}
    paths::Union{Nothing, Dict{String, String}}

    function EarthBoxPathsState()
        key_names = EarthBoxPathNames()
        required_paths = get_required_paths(key_names)
        new(key_names, required_paths, nothing)
    end
end

function get_required_paths(key_names::EarthBoxPathNames)
    return [
        key_names.model_input_file,
        key_names.materials_input_file,
        key_names.materials_library_file,
        key_names.output_dir
    ]
end

function set_paths!(
    eb_paths_state::EarthBoxPathsState, 
    input_paths::Union{Nothing, Dict{String, String}}=nothing
)
    if isnothing(input_paths)
        input_paths = Dict{String, String}()
    end

    for name in eb_paths_state.required_paths
        if !haskey(input_paths, name)
            if name != "output_dir"
                input_paths = update_input_paths_dict(input_paths, name, nothing)
            else
                input_paths = update_input_paths_dict(input_paths, name, name)
            end
        end
    end

    input_paths["active_model_dir"] = pwd()
    eb_paths_state.paths = input_paths
end

function update_input_paths_dict(
    input_paths::Dict{String, String},
    path_name::String,
    updated_value::Union{String, Nothing}
)::Dict{String, String}
    if !haskey(input_paths, path_name)
        print_info("A $(path_name) was not defined in the paths dictionary. The user will have to provide inputs through the earthbox API.")
        input_paths[path_name] = string(updated_value)
    end
    return input_paths
end

function check_output_dir(paths::EarthBoxPathsState)
    if !isdir(paths.paths["output_dir"])
        print_info("Output directory not found. A directory will be created for you.")
        mkpath(paths.paths["output_dir"])
    end
    print_info("Output directory: $(paths.paths["output_dir"])")
end

function print_input_path_info(eb_paths_state::EarthBoxPathsState)
    print_info("Input path information at EarthBox initialization:", level=1)
    print_info("Active model directory: $(eb_paths_state.paths["active_model_dir"])", level=2)
    print_info("Model input file: $(eb_paths_state.paths["model_input_file"])", level=2)
    print_info("Materials input file: $(eb_paths_state.paths["materials_input_file"])", level=2)
    print_info("Material library file: $(eb_paths_state.paths["materials_library_file"])", level=2)
    print_info("Output directory: $(eb_paths_state.paths["output_dir"])", level=2)
end

function create_run_paths_dict(
    eb_paths::EarthBoxPathsState, 
    base_path::String
)::Union{Nothing, Dict{String, String}}
    if !isnothing(eb_paths.paths)
        run_paths = Dict{String, String}()
        run_paths["model_file"] = eb_paths.paths["model_input_file"]
        run_paths["materials_file"] = eb_paths.paths["materials_input_file"]
        run_paths["materials_library_file"] = eb_paths.paths["materials_library_file"]
        run_paths["output_dir"] = eb_paths.paths["output_dir"]
        run_paths["src_dir"] = base_path
        run_paths["plots_dir"] = joinpath(eb_paths.paths["output_dir"], "plots")
        return run_paths
    end
    return nothing
end

end # module 