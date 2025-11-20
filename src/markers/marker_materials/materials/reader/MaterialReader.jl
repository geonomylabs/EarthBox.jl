module MaterialReader

import YAML
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.EarthBoxDtypes: MaterialDictType
import EarthBox.EarthBoxDtypes: RawMaterialParameters
import EarthBox.EarthBoxDtypes: RawMaterialsInputType
import EarthBox.EarthBoxDtypes: InputMaterialParamDictType
import EarthBox.ParameterRegistry: get_master_material_parameters

function build_materials_dict(
    materials_library_dict::MaterialsDictType,
    materials_input_dict::MaterialsDictType
)::Tuple{MaterialsDictType, Int64}
    materials_dict = load_rock_properties_from_library(
        materials_library_dict, materials_input_dict)
    nmats = length(materials_dict)
    return materials_dict, nmats
end

function load_rock_properties_from_library(
    materials_library::MaterialsDictType,
    materials_input_dict::MaterialsDictType
)::MaterialsDictType
    materials_dict = deepcopy(materials_input_dict)
    for (matid, parameters) in materials_dict
        material_name = parameters["mat_name"]
        parameters_library = materials_library[material_name]
        for (param_name, param_value) in parameters_library
            materials_dict[matid][param_name] = param_value
        end
    end
    return materials_dict
end

function read_materials_library(
    material_library_file_path::String
)::Tuple{MaterialsDictType, Int64}
    file_type = "library"
    (materials_library_dict, nmats_lib) = read_material_file(material_library_file_path, file_type)
    return materials_library_dict, nmats_lib
end

function read_materials_input(
    material_model_file_path::String
)::Tuple{MaterialsDictType, Int64}
    file_type = "model"
    (materials_dict, nmats) = read_material_file(material_model_file_path, file_type)
    return materials_dict, nmats
end

function read_material_file(
    material_file_path::String,
    file_type::String
)::Tuple{MaterialsDictType, Int64}
    raw_material_file_data = read_yaml_file(file_type, material_file_path)
    _, master_parameters = get_master_material_parameters()
    materials_dict = MaterialsDictType()
    for (material_name, raw_parameters_dict) in raw_material_file_data
        if isa(material_name, Int64)
            material_name = convert(Int16, material_name)
        end
        if !isa_header(material_name)
            material_parameters_dict = create_parameters_dict(
                file_type, material_name, raw_parameters_dict, master_parameters)
            # Convert material_name to the correct type for MaterialsDictType
            material_key = convert(Union{String, Int16}, material_name)
            materials_dict[material_key] = material_parameters_dict
        end
    end
    return materials_dict, length(keys(materials_dict))
end

""" Read a YAML file and return a dictionary.

# Material Library Yaml File Format

title: "Title of the material collection"
description: "Description of the material collection"
nmats: Number of materials in the collection
Material_Collection_Name_1:
    Material_Name_1:
        parameter_name_1: [parameter_value_1, parameter_units_1, parameter_description_1]
        parameter_name_2: [parameter_value_2, parameter_units_2, parameter_description_2]
    Material_Name_2:
        parameter_name_1: [parameter_value_1, parameter_units_1, parameter_description_1]
        parameter_name_2: [parameter_value_2, parameter_units_2, parameter_description_2]
Material_Collection_Name_2:
    Material_Name_3:
        parameter_name_1: [parameter_value_1, parameter_units_1, parameter_description_1]
        parameter_name_2: [parameter_value_2, parameter_units_2, parameter_description_2]
    Material_Name_4:
        parameter_name_1: [parameter_value_1, parameter_units_1, parameter_description_1]
        parameter_name_2: [parameter_value_2, parameter_units_2, parameter_description_2]

Note that material collection names are filtered out of the final dictionary.
Therefore, material names must be unique across material collections.

# Material Input Yaml File Format

name: Material name
material_model_id: Material model ID
Material_ID_0:
    parameter_name_1: [parameter_value_1, parameter_units_1, parameter_description_1]
    parameter_name_2: [parameter_value_2, parameter_units_2, parameter_description_2]
Material_ID_1:
    parameter_name_1: [parameter_value_1, parameter_units_1, parameter_description_1]
    parameter_name_2: [parameter_value_2, parameter_units_2, parameter_description_2]
Material_ID_2:
    parameter_name_1: [parameter_value_1, parameter_units_1, parameter_description_1]
    parameter_name_2: [parameter_value_2, parameter_units_2, parameter_description_2]
Material_ID_2:
    parameter_name_1: [parameter_value_1, parameter_units_1, parameter_description_1]
    parameter_name_2: [parameter_value_2, parameter_units_2, parameter_description_2]

# Returns
- data_dict::Dict{Any, Any}: Dictionary where keys are header or material names and
    values are header values or material data, respectively.
"""
function read_yaml_file(
    file_type::String,
    yaml_file_path::String
)::Dict{Any, Any}
    if isfile(yaml_file_path)
        data_dict = YAML.load_file(yaml_file_path)
        if file_type == "library"
            data_dict = remove_collection_names_from_library_dictionary(data_dict)
        end
        return data_dict
    end
    throw(ArgumentError("File $yaml_file_path does not exist."))
end

function remove_collection_names_from_library_dictionary(
    data_dict::Dict{Any, Any}
)::Dict{Any, Any}
    data_dict_new = Dict{Any, Any}()
    name_list = String[]
    for (collection_name, collection_data) in data_dict
        if isa(collection_name, Int64)
            collection_name = convert(Int16, collection_name)
        end
        if !isa_header(collection_name)
            for (material_name, material_data) in collection_data
                if material_name ∉ name_list
                    push!(name_list, material_name)
                else
                    throw(ArgumentError(get_duplicate_error_msg(material_name)))
                end
                data_dict_new[material_name] = material_data
            end
        else
            data_dict_new[collection_name] = collection_data
        end
    end
    return data_dict_new
end

function isa_header(section::Union{String, Int16})::Bool
    return section in ("title", "name", "material_model_id", "description", "nmats")
end

function get_duplicate_error_msg(material_name::String)::String
    return (
        "Found duplicate material name in library: $material_name. "
        *"All material names must be unique across material collections "
        *"in a material library file."
    )
end

function create_parameters_dict(
    file_type::String,
    material_identifier::Union{String, Int16},
    raw_parameters_dict::Dict{Any, Any},
    master_parameters::InputMaterialParamDictType
)::MaterialDictType
    parameters_dict = MaterialDictType()
    parameters_dict[get_parameter_name_for_identifier(file_type)] = material_identifier
    for (parameter_name, parameter_data_vector) in raw_parameters_dict
        check_units(parameter_name, parameter_data_vector[2], master_parameters)
        parameters_dict[parameter_name] = parameter_data_vector[1]
    end
    if file_type == "library"
        add_missing_parameters_from_library(parameters_dict, master_parameters)
    end
    parameters_dict = convert_strings_to_master_type(parameters_dict, master_parameters)
    return parameters_dict
end

""" Get the parameter name for the material identifier.

If the file type is "model", the material identifier is "matid". Otherwise, the
file type is "library" and the material identifier is "name".

"""
function get_parameter_name_for_identifier(file_type::String)::String
    return file_type == "model" ? "matid" : "name"
end

function add_missing_parameters_from_library(
    param_dict::MaterialDictType,
    master_parameters::InputMaterialParamDictType
)::Nothing
    current_param_names = collect(keys(param_dict))
    for (param_name_master, param_info_list) in master_parameters
        if applicable_to_library(param_name_master)
            if param_name_master ∉ current_param_names
                param_dict[param_name_master] = param_info_list[2]
            end
        end
    end
    return nothing
end

function applicable_to_library(param_name_master::String)::Bool
    names_from_model_material_file = [
        "mat_name",
        "mat_type",
        "mat_domain",
        "red_fraction",
        "green_fraction",
        "blue_fraction"
    ]
    header_names = ["nmats", "description", "material_model_id", "title"]
    obsolete_names = ["ncolors", "matid", "itype_mat"]
    return !(param_name_master in names_from_model_material_file ||
             param_name_master in header_names ||
             param_name_master in obsolete_names)
end

function convert_strings_to_master_type(
    parameters_dict::MaterialDictType,
    master_parameters::InputMaterialParamDictType
)::MaterialDictType
    master_names = collect(keys(master_parameters))
    parameters_dict_new = MaterialDictType()
    for (parameter_name, parameter_value) in parameters_dict
        check_name(parameter_name, master_names)
        parameter_value = convert_parameter_to_data_type(
            parameter_value, parameter_name, master_parameters)
        parameters_dict_new[parameter_name] = parameter_value
    end
    return parameters_dict_new
end

function check_name(param_name::String, master_names::Vector{String})::Nothing
    if param_name ∉ master_names
        throw(ArgumentError("Name $param_name is invalid."))
    end
    return nothing
end

function convert_parameter_to_data_type(
    parameter_value::Union{String, Int16, Int64, Float64},
    parameter_name::String,
    master_parameters::InputMaterialParamDictType
)::Union{String, Int16, Int64, Float64}
    conversion_func = master_parameters[parameter_name][1]
    if isa(parameter_value, String)
        return conversion_func(parameter_value)
    else
        return parameter_value
    end
end

function check_units(
    param_name::String,
    param_units::String,
    master_parameters::InputMaterialParamDictType
)::Nothing
    master_units = master_parameters[param_name][3]
    if param_units != master_units
        throw(ArgumentError(
            "Units $param_units do not match master units $master_units for parameter $param_name."
        ))
    end
    return nothing
end

end # module 