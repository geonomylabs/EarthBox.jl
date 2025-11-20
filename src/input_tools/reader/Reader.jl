module Reader

import YAML
import EarthBox.EarthBoxDtypes: ParametersDictType

struct InputData
    dfloats::Dict{String, Float64}
    dints::Dict{String, Int64}
    dstrs::Dict{String, String}
end

function get_parameters_input_dict(
    input_file_path::Union{String, Nothing}
)::Union{Dict{String, Vector{Any}}, Nothing}
    parameters_input_dict = nothing
    if input_file_path !== nothing && isfile(input_file_path)
        model_input_dict = read_yaml_file(input_file_path)
        parameters_input_dict = make_parameters_dict(model_input_dict)
    end
    return parameters_input_dict
end

function get_input_data(
    input_file_path::Union{String, Nothing}, 
    initialization_params::Union{Dict{String, Union{Float64, Int64}}, Nothing}
)
    input_data_obj = build_input_data_structure()
    initialize_input_data_structure(input_data_obj)

    if input_file_path !== nothing && isfile(input_file_path)
        load_input_file_parameters_into_input_data_structure(
            input_file_path, input_data_obj)
    end
    if initialization_params !== nothing
        load_parameter_dict(initialization_params, input_data_obj)
    end
    return input_data_obj
end


function initialize_input_data_structure(input_data_obj::InputData)
    _, input_param_dict = get_input_parameter_name_list()
    for (param_name, info_list) in input_param_dict
        param_value = info_list[2]
        load_inputdata_value(input_data_obj, param_value, param_name)
    end
end

function read_yaml_file(yaml_file_path::String)
    data_dict = Dict()
    if yaml_file_path !== nothing
        data_dict = YAML.load_file(yaml_file_path)
    end
    return data_dict
end

function load_input_file_parameters_into_input_data_structure(
    input_file_path::String, input_data_obj::InputData
)
    model_input_dict = read_yaml_file(input_file_path)
    parameters_dict = make_parameters_dict(model_input_dict)
    load_parameter_dict(parameters_dict, input_data_obj)
end

"""
    make_parameters_dict(model_input_dict::Dict)

Read inputs from yaml file and remove section keys.

# Input
- `model_input_dict` : Dict
    Input dictionary read from yaml file that contains section keys and
    input parameter lists.

# Returns
- `parameters_dict` : Dict{String, Vector{Any}}
    Dictionary with input parameter name keys and elements equal to
    input parameter lists where each list includes parameter value, units
    and description.

"""
function make_parameters_dict(model_input_dict::Dict)::Dict{String, Vector{Any}}
    parameters_dict = Dict{String, Vector{Any}}()
    for (section_name, sect_element) in model_input_dict
        if sect_element isa Dict
            for (param_name, param_list) in sect_element
                parameters_dict[param_name] = param_list
            end
        elseif sect_element isa Vector
            parameters_dict[section_name] = sect_element
        end
    end
    return parameters_dict
end

"""
    load_parameter_dict(parameters_dict::Dict{String, Vector{Any}}, input_data_obj::InputData)

Load parameter dictionary into input data object.
"""
function load_parameter_dict(
    parameters_dict::Dict{String, Vector{Any}}, 
    input_data_obj::InputData
)
    input_param_names, input_param_dict = get_input_parameter_name_list()
    for (param_name, param_list) in parameters_dict
        if param_name ∉ input_param_names
            throw(ArgumentError("Name $param_name in yaml file is invalid."))
        end
        type_conversion_func = input_param_dict[param_name][1]
        value = type_conversion_func(param_list[1])
        load_inputdata_value(input_data_obj, value, param_name)
    end
end

"""
    load_inputdata_value(data_obj::InputData, param_value::Union{Float64, Int64, String}, param_name::String)

Load input data into data object dictionaries.
"""
function load_inputdata_value(
    data_obj::InputData, 
    param_value::Union{Float64, Int64, String},
    param_name::String
)
    if param_value isa Float64
        data_obj.dfloats[param_name] = param_value
    elseif param_value isa Int64
        data_obj.dints[param_name] = param_value
    elseif param_value isa String
        data_obj.dstrs[param_name] = param_value
    end
end

"""
    get_parameter_info(elements::Vector{String}, input_param_names::Vector{String}, iline::Int)

Get parameter information from file line elements.
"""
function get_parameter_info(
    elements::Vector{String}, 
    input_param_names::Vector{String}, 
    iline::Int
)
    param_value_str = split(elements[1])
    param_name = split(elements[2])
    
    if !isempty(param_name)
        param_value_str = replace(param_value_str[1], " " => "")
        param_name = replace(param_name[1], " " => "")
        if param_name ∉ input_param_names
            throw(ArgumentError(
                "Attribute name $param_name is invalid on line $(iline+1) of "
                *"input file.")
                )
        end
    else
        throw(ArgumentError(
            "Attribute name $param_name is invalid on line $(iline+1) of input "
            *"file. Length of variable name is zero.")
            )
    end
    
    return param_value_str, param_name
end

end # module 