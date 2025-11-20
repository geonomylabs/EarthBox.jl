module InternalVelocityZoneGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "internal_velocity_zone"

const PDATA = get_eb_parameters()

"""
    InternalVelocityZone <: AbstractParameterGroup

Parameter group for internal velocity zone geometry parameters.

# Fields
- `xindex_vx_internal::`[`ParameterInt`](@ref): $(PDATA.xindex_vx_internal.description)
- `yindex_min_vx_internal::`[`ParameterInt`](@ref): $(PDATA.yindex_min_vx_internal.description)
- `yindex_max_vx_internal::`[`ParameterInt`](@ref): $(PDATA.yindex_max_vx_internal.description)
- `xindex_vy_internal::`[`ParameterInt`](@ref): $(PDATA.xindex_vy_internal.description)
- `yindex_min_vy_internal::`[`ParameterInt`](@ref): $(PDATA.yindex_min_vy_internal.description)
- `yindex_max_vy_internal::`[`ParameterInt`](@ref): $(PDATA.yindex_max_vy_internal.description)
- `x_vx_internal::`[`ParameterFloat`](@ref): $(PDATA.x_vx_internal.description)
- `y_min_vx_internal::`[`ParameterFloat`](@ref): $(PDATA.y_min_vx_internal.description)
- `y_max_vx_internal::`[`ParameterFloat`](@ref): $(PDATA.y_max_vx_internal.description)
- `x_vy_internal::`[`ParameterFloat`](@ref): $(PDATA.x_vy_internal.description)
- `y_min_vy_internal::`[`ParameterFloat`](@ref): $(PDATA.y_min_vy_internal.description)
- `y_max_vy_internal::`[`ParameterFloat`](@ref): $(PDATA.y_max_vy_internal.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
  parameter objects

# Nested Dot Access
- `xindex_vx_internal = $(ROOT_NAME).$(GRP_NAME).xindex_vx_internal.value`
- `yindex_min_vx_internal = $(ROOT_NAME).$(GRP_NAME).yindex_min_vx_internal.value`
- `yindex_max_vx_internal = $(ROOT_NAME).$(GRP_NAME).yindex_max_vx_internal.value`
- `xindex_vy_internal = $(ROOT_NAME).$(GRP_NAME).xindex_vy_internal.value`
- `yindex_min_vy_internal = $(ROOT_NAME).$(GRP_NAME).yindex_min_vy_internal.value`
- `yindex_max_vy_internal = $(ROOT_NAME).$(GRP_NAME).yindex_max_vy_internal.value`
- `x_vx_internal = $(ROOT_NAME).$(GRP_NAME).x_vx_internal.value`
- `y_min_vx_internal = $(ROOT_NAME).$(GRP_NAME).y_min_vx_internal.value`
- `y_max_vx_internal = $(ROOT_NAME).$(GRP_NAME).y_max_vx_internal.value`
- `x_vy_internal = $(ROOT_NAME).$(GRP_NAME).x_vy_internal.value`
- `y_min_vy_internal = $(ROOT_NAME).$(GRP_NAME).y_min_vy_internal.value`
- `y_max_vy_internal = $(ROOT_NAME).$(GRP_NAME).y_max_vy_internal.value`

# Constructor
    InternalVelocityZone()

Create a new InternalVelocityZone parameter group with default values.

# Returns
- `InternalVelocityZone`: New InternalVelocityZone parameter group with initialized default 
  values

"""
mutable struct InternalVelocityZone <: AbstractParameterGroup
    xindex_vx_internal::ParameterInt
    yindex_min_vx_internal::ParameterInt
    yindex_max_vx_internal::ParameterInt
    xindex_vy_internal::ParameterInt
    yindex_min_vy_internal::ParameterInt
    yindex_max_vy_internal::ParameterInt
    x_vx_internal::ParameterFloat
    y_min_vx_internal::ParameterFloat
    y_max_vx_internal::ParameterFloat
    x_vy_internal::ParameterFloat
    y_min_vy_internal::ParameterFloat
    y_max_vy_internal::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function InternalVelocityZone()::InternalVelocityZone
    pdata = get_eb_parameters()
    data = InternalVelocityZone(
        pdata.xindex_vx_internal,
        pdata.yindex_min_vx_internal,
        pdata.yindex_max_vx_internal,
        pdata.xindex_vy_internal,
        pdata.yindex_min_vy_internal,
        pdata.yindex_max_vy_internal,
        pdata.x_vx_internal,
        pdata.y_min_vx_internal,
        pdata.y_max_vx_internal,
        pdata.x_vy_internal,
        pdata.y_min_vy_internal,
        pdata.y_max_vy_internal,
        Union{ParameterFloat}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 