module SetParameters

import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.Parameters: set_parameter_float!, set_parameter_int!

function check_and_update_parameters(
    obj_list::Vector{Union{ParameterFloat, ParameterInt}},
    parameters::Dict{String, <:Union{Float64, Int64, String}}
)::Nothing
    parameter_object_dict = get_parameter_object_dict(obj_list)
    parameter_object_names = collect(keys(parameter_object_dict))
    already_set_names = String[]
    for (name, value) in parameters
        if name ∉ parameter_object_names
            invalid_parameter_name(name)
        end
        if name in already_set_names
            duplicate_parameter_name(name)
        end
        set_parameter(parameter_object_dict[name], value)
        push!(already_set_names, name)
    end
    return nothing
end

function get_parameter_object_dict(
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
)::Dict{String, Union{ParameterFloat, ParameterInt}}
    master_dict = Dict{String, Union{ParameterFloat, ParameterInt}}()
    for obj in obj_list
        master_dict[obj.name] = obj
    end
    return master_dict
end

function set_parameter(
    param_obj::Union{ParameterFloat, ParameterInt},
    value::Union{Float64, Int64}
)::Nothing
    if param_obj isa ParameterFloat && value isa Float64
        set_parameter_float!(param_obj, value)
    elseif param_obj isa ParameterInt && value isa Int64
        set_parameter_int!(param_obj, value)
    end
    return nothing
end

function check_and_get_parameter_value(
    obj_list::Vector{Union{ParameterFloat, ParameterInt}},
    parameter_name::String
)::Union{Float64, Int64, String}
    master_dict = Dict{String, Union{ParameterFloat, ParameterInt}}()
    for obj in obj_list
        master_dict[obj.name] = obj
    end
    if haskey(master_dict, parameter_name)
        return master_dict[parameter_name].value
    else
        invalid_parameter_name(parameter_name)
    end
end

function invalid_parameter_name(name::String)::Nothing
    output_str = (
        "Input parameter name $(name) is not valid. Check your inputs.")
    throw(ArgumentError(output_str))
    return nothing
end

function duplicate_parameter_name(name::String)::Nothing
    output_str = (
        "Input parameter name $(name) was encountered twice. "
        *"Check your inputs."
        )
    throw(ArgumentError(output_str))
    return nothing
end

end # module
