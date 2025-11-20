module PlumeGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "plume"

const PDATA = get_eb_parameters()

"""
    Plume <: AbstractParameterGroup

Parameter group for mantle plume geometry parameters.

# Fields
- `plume_radius::`[`ParameterFloat`](@ref): $(PDATA.plume_radius.description)
- `plume_center_x::`[`ParameterFloat`](@ref): $(PDATA.plume_center_x.description)
- `plume_center_y::`[`ParameterFloat`](@ref): $(PDATA.plume_center_y.description)
- `plume_head_thick::`[`ParameterFloat`](@ref): $(PDATA.plume_head_thick.description)
- `delta_temperature_plume::`[`ParameterFloat`](@ref): $(PDATA.delta_temperature_plume.description)
- `iuse_plume::`[`ParameterInt`](@ref): $(PDATA.iuse_plume.description)
- `plume_start_time::`[`ParameterFloat`](@ref): $(PDATA.plume_start_time.description)
- `iplume_state::`[`ParameterInt`](@ref): $(PDATA.iplume_state.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
    parameter objects

# Nested Dot Access
- `plume_radius = $(ROOT_NAME).$(GRP_NAME).plume_radius.value`
- `plume_center_x = $(ROOT_NAME).$(GRP_NAME).plume_center_x.value`
- `plume_center_y = $(ROOT_NAME).$(GRP_NAME).plume_center_y.value`
- `plume_head_thick = $(ROOT_NAME).$(GRP_NAME).plume_head_thick.value`
- `delta_temperature_plume = $(ROOT_NAME).$(GRP_NAME).delta_temperature_plume.value`
- `iuse_plume = $(ROOT_NAME).$(GRP_NAME).iuse_plume.value`
- `plume_start_time = $(ROOT_NAME).$(GRP_NAME).plume_start_time.value`
- `iplume_state = $(ROOT_NAME).$(GRP_NAME).iplume_state.value`

# Constructor
    Plume()

Create a new Plume parameter group with default values.

# Returns
- `Plume`: New Plume parameter group with initialized default values

"""
mutable struct Plume <: AbstractParameterGroup
    plume_radius::ParameterFloat
    plume_center_x::ParameterFloat
    plume_center_y::ParameterFloat
    plume_head_thick::ParameterFloat
    delta_temperature_plume::ParameterFloat
    iuse_plume::ParameterInt
    plume_start_time::ParameterFloat
    iplume_state::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Plume()::Plume
    pdata = get_eb_parameters()
    data = Plume(
        pdata.plume_radius,
        pdata.plume_center_x,
        pdata.plume_center_y,
        pdata.plume_head_thick,
        pdata.delta_temperature_plume,
        pdata.iuse_plume,
        pdata.plume_start_time,
        pdata.iplume_state,
        Union{ParameterFloat, ParameterInt}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 