module LithoStrongZonesGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "litho_strong_zones"

const PDATA = get_eb_parameters()

"""
    LithoStrongZones <: AbstractParameterGroup

Parameter group for lithospheric strong zones geometry parameters.

# Fields
- `x_left_strong::`[`ParameterFloat`](@ref): $(PDATA.x_left_strong.description)
- `x_right_strong::`[`ParameterFloat`](@ref): $(PDATA.x_right_strong.description)
- `iuse_strong_zones::`[`ParameterInt`](@ref): $(PDATA.iuse_strong_zones.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter
    objects

# Nested Dot Access
- `x_left_strong = $(ROOT_NAME).$(GRP_NAME).x_left_strong.value`
- `x_right_strong = $(ROOT_NAME).$(GRP_NAME).x_right_strong.value`
- `iuse_strong_zones = $(ROOT_NAME).$(GRP_NAME).iuse_strong_zones.value`

# Constructor
    LithoStrongZones()

Create a new LithoStrongZones parameter group with default values.

# Returns
- `LithoStrongZones`: New LithoStrongZones parameter group with initialized default values

"""
mutable struct LithoStrongZones <: AbstractParameterGroup
    x_left_strong::ParameterFloat
    x_right_strong::ParameterFloat
    iuse_strong_zones::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function LithoStrongZones()::LithoStrongZones
    pdata = get_eb_parameters()
    data = LithoStrongZones(
        pdata.x_left_strong,
        pdata.x_right_strong,
        pdata.iuse_strong_zones,
        Union{ParameterFloat, ParameterInt}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 