module ParameterGroupTools

import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.SetParameters: check_and_update_parameters
import EarthBox.SetParameters: check_and_get_parameter_value
import EarthBox.EarthBoxDtypes: InputDictType
import EarthBox.EarthBoxDtypes: AbstractParameterGroup
import EarthBox.EarthBoxDtypes: AbstractParameterCollection

function set_group_parameters!(
    parameter_group::AbstractParameterGroup,
    parameters::Dict{String, <:Union{Float64, Int64, String}}
)::Nothing
    check_and_update_parameters(parameter_group.obj_list, parameters)
    return nothing
end

function get_parameter_value(
    parameter_group::AbstractParameterGroup,
    parameter_name::String
)::Union{Float64, Int, String}
    return check_and_get_parameter_value(parameter_group.obj_list, parameter_name)
end

function print(parameter_group::AbstractParameterGroup)
    for obj in parameter_group.obj_list
        print_parameter(obj)
    end
end

function get_numerical_parameter_object_list(
    data::Union{AbstractParameterGroup, AbstractParameterCollection}
)::Vector{Union{ParameterFloat, ParameterInt}}
    obj_list = []
    for field in fieldnames(typeof(data))
        val = getfield(data, field)
        if val isa ParameterFloat || val isa ParameterInt
            push!(obj_list, val)
        end
    end
    return obj_list
end

end # module