module GetPaths

import EarthBox.GetArgs: get_model_output_path_from_args, get_storage_path_from_args

"""
    get_storage_path(
        model_case_name::String,
        base_path::String
    )::String

Get the storage path from the command line argument `storage_path=...`. If the 
storage path is not provided, then the storage path is generated from the model 
case name and base path.

# Arguments
- `model_case_name::String`:
    - Name of the model case
- `root_path::String`:
    - Root path of the storage directory

# Returns
- `storage_path::String`:
    - Path to the storage directory
"""
function get_storage_path(
    model_case_name::String,
    root_path::Union{String, Nothing}=nothing
)::String
    storage_path = get_storage_path_from_args()
    if isnothing(storage_path)
        if isnothing(root_path)
            throw(ArgumentError(
                "storage_path is not provided and root_path is not provided. "
                *"A path to a storage directory cannot be defined.")
                )
        end
        storage_path = joinpath(root_path, "$(model_case_name)_output")
    end
    if !isdir(storage_path)
        mkpath(storage_path)
    end
    return storage_path
end

"""
    get_model_output_path(
        model_case_name::String,
        base_path::String
    )::String

Get the model output path or storage path from the command line argument 
`model_output_path=...`. If the model output path is not provided, then the 
model output path is generated from the model case name and base path.

# Arguments
- `model_case_name::String`:
    - Name of the model case
- `root_path::String`:
    - Root path of the model output directory

# Returns
- `model_output_path::String`:
    - Path to the model output directory
"""
function get_model_output_path(
    model_case_name::String,
    root_path::Union{String, Nothing}=nothing
)::String
    model_output_path = get_model_output_path_from_args()
    if isnothing(model_output_path)
        if isnothing(root_path)
            throw(ArgumentError(
                "model_output_path is not provided and root_path is not provided. "
                *"A path to a model output directory cannot be defined.")
                )
        end
        model_output_path = joinpath(root_path, "$(model_case_name)_output")
    end
    if !isdir(model_output_path)
        mkpath(model_output_path)
    end
    return model_output_path
end

function get_path_to_materials_library()
    parent = parentmodule(GetPaths)
    src_dir = dirname(pathof(parent))
    local_library_path = get_local_material_library_path()
    return joinpath(src_dir, local_library_path)
end

function get_earthbox_src_path()
    parent = parentmodule(GetPaths)
    return dirname(pathof(parent))
end

function get_local_material_library_path()
    return "materials_library/libraries"
end

function get_output_path(
    drive_root_path::String,
    drive_base_name::String,
    drive_number_id::Int,
    local_group_path_on_drive::String,
    model_base_name::String,
    model_case_name::String
)
    experiment_group_dir_path = get_experiment_group_dir_path(
        drive_root_path, drive_base_name, drive_number_id,
        local_group_path_on_drive
    )
    model_name = model_base_name * "_" * model_case_name
    output_dir_name = model_name * "_output"
    output_dir_path = joinpath(experiment_group_dir_path, output_dir_name)
    return output_dir_path
end

function get_experiment_group_dir_path(
    drive_root_path::String,
    drive_base_name::String,
    drive_number_id::Int,
    local_group_path_on_drive::String
)
    drive_path = get_drive_path(drive_root_path, drive_base_name, drive_number_id)
    experiment_group_dir_path = joinpath(drive_path, local_group_path_on_drive)
    if !isdir(experiment_group_dir_path)
        mkpath(experiment_group_dir_path)
    end
    return experiment_group_dir_path
end

function get_drive_path(
    drive_root_path::String,
    drive_base_name::String,
    drive_number_id::Int
)
    drive_name = get_drive_name(drive_base_name, drive_number_id)
    drive_path = joinpath(drive_root_path, drive_name)
    return drive_path
end

function get_drive_name(
    drive_base_name::String,
    drive_number::Int
)
    return drive_base_name * string(drive_number)
end

end # module GetPaths 