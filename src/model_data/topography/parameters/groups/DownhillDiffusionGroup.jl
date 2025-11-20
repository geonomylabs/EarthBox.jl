module DownhillDiffusionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.topography.parameters"
const GRP_NAME = "downhill_diffusion"

const PDATA = get_eb_parameters()

"""
    DownhillDiffusion <: AbstractParameterGroup

Parameter group for downhill diffusion parameters.

# Fields
- `downhill_diff_elev_max::`[`ParameterFloat`](@ref): $(PDATA.downhill_diff_elev_max.description)
- `transport_length::`[`ParameterFloat`](@ref): $(PDATA.transport_length.description)
- `topo_diff_coef::`[`ParameterFloat`](@ref): $(PDATA.topo_diff_coef.description)
- `subaerial_slope_diffusivity::`[`ParameterFloat`](@ref): $(PDATA.subaerial_slope_diffusivity.description)
- `precipitation_rate::`[`ParameterFloat`](@ref): $(PDATA.precipitation_rate.description)
- `subaerial_transport_coefficient::`[`ParameterFloat`](@ref): $(PDATA.subaerial_transport_coefficient.description)
- `submarine_slope_diffusivity::`[`ParameterFloat`](@ref): $(PDATA.submarine_slope_diffusivity.description)
- `submarine_diffusion_decay_depth::`[`ParameterFloat`](@ref): $(PDATA.submarine_diffusion_decay_depth.description)
- `transport_timestep::`[`ParameterFloat`](@ref): $(PDATA.transport_timestep.description)
- `number_of_transport_timesteps_per_model_timestep::`[`ParameterInt`](@ref): $(PDATA.number_of_transport_timesteps_per_model_timestep.description)
- `iuse_compaction_correction::`[`ParameterInt`](@ref): $(PDATA.iuse_compaction_correction.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `downhill_diff_elev_max = $(ROOT_NAME).$(GRP_NAME).downhill_diff_elev_max.value`
- `transport_length = $(ROOT_NAME).$(GRP_NAME).transport_length.value`
- `topo_diff_coef = $(ROOT_NAME).$(GRP_NAME).topo_diff_coef.value`
- `subaerial_slope_diffusivity = $(ROOT_NAME).$(GRP_NAME).subaerial_slope_diffusivity.value`
- `precipitation_rate = $(ROOT_NAME).$(GRP_NAME).precipitation_rate.value`
- `subaerial_transport_coefficient = $(ROOT_NAME).$(GRP_NAME).subaerial_transport_coefficient.value`
- `submarine_slope_diffusivity = $(ROOT_NAME).$(GRP_NAME).submarine_slope_diffusivity.value`
- `submarine_diffusion_decay_depth = $(ROOT_NAME).$(GRP_NAME).submarine_diffusion_decay_depth.value`
- `transport_timestep = $(ROOT_NAME).$(GRP_NAME).transport_timestep.value`
- `iuse_compaction_correction = $(ROOT_NAME).$(GRP_NAME).iuse_compaction_correction.value`

# Constructor
    DownhillDiffusion()

# Returns
- `DownhillDiffusion`: New DownhillDiffusion parameter group with initialized values

"""
mutable struct DownhillDiffusion <: AbstractParameterGroup
    downhill_diff_elev_max::ParameterFloat
    transport_length::ParameterFloat
    topo_diff_coef::ParameterFloat
    subaerial_slope_diffusivity::ParameterFloat
    precipitation_rate::ParameterFloat
    subaerial_transport_coefficient::ParameterFloat
    submarine_slope_diffusivity::ParameterFloat
    submarine_diffusion_decay_depth::ParameterFloat
    transport_timestep::ParameterFloat
    number_of_transport_timesteps_per_model_timestep::ParameterInt
    iuse_compaction_correction::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function DownhillDiffusion()::DownhillDiffusion
    pdata = get_eb_parameters()
    data = DownhillDiffusion(
        pdata.downhill_diff_elev_max,
        pdata.transport_length,
        pdata.topo_diff_coef,
        pdata.subaerial_slope_diffusivity,
        pdata.precipitation_rate,
        pdata.subaerial_transport_coefficient,
        pdata.submarine_slope_diffusivity,
        pdata.submarine_diffusion_decay_depth,
        pdata.transport_timestep,
        pdata.number_of_transport_timesteps_per_model_timestep,
        pdata.iuse_compaction_correction,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
