module TempChangeLimitGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.heat_equation.parameters"
const GRP_NAME = "temp_change_limit"

const PDATA = get_eb_parameters()

"""
    TempChangeLimit <: AbstractParameterGroup

Parameter group for temperature change limits.

# Fields
- `max_temp_change::`[`ParameterFloat`](@ref): $(PDATA.max_temp_change.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `max_temp_change = $(ROOT_NAME).$(GRP_NAME).max_temp_change.value`

# Constructor
    TempChangeLimit()

# Returns
- `TempChangeLimit`: New TempChangeLimit parameter group with initialized values

"""
mutable struct TempChangeLimit <: AbstractParameterGroup
    max_temp_change::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function TempChangeLimit()::TempChangeLimit
    pdata = get_eb_parameters()
    data = TempChangeLimit(
        pdata.max_temp_change,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
