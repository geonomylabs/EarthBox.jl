module HydrothermalGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "hydrothermal"

const PDATA = get_eb_parameters()

"""
    Hydrothermal <: AbstractParameterGroup

Parameter group for hydrothermal circulation properties.

# Fields
- `iuse_hydrothermal::`[`ParameterInt`](@ref): $(PDATA.iuse_hydrothermal.description)
- `iuse_melt_lens::`[`ParameterInt`](@ref): $(PDATA.iuse_melt_lens.description)
- `hydrothermal_smoothing_factor::`[`ParameterFloat`](@ref): $(PDATA.hydrothermal_smoothing_factor.description)
- `hydrothermal_nusselt_number::`[`ParameterFloat`](@ref): $(PDATA.hydrothermal_nusselt_number.description)
- `hydrothermal_max_temperature::`[`ParameterFloat`](@ref): $(PDATA.hydrothermal_max_temperature.description)
- `hydrothermal_max_depth::`[`ParameterFloat`](@ref): $(PDATA.hydrothermal_max_depth.description)
- `iuse_plastic_strain_rate_for_hydrothermal::`[`ParameterInt`](@ref): $(PDATA.iuse_plastic_strain_rate_for_hydrothermal.description)
- `hydrothermal_decay_length::`[`ParameterFloat`](@ref): $(PDATA.hydrothermal_decay_length.description)
- `hydrothermal_buffer_distance::`[`ParameterFloat`](@ref): $(PDATA.hydrothermal_buffer_distance.description)
- `hydrothermal_plastic_strain_rate_reference::`[`ParameterFloat`](@ref): $(PDATA.hydrothermal_plastic_strain_rate_reference.description)
- `iuse_plastic_strain_for_hydrothermal::`[`ParameterInt`](@ref): $(PDATA.iuse_plastic_strain_for_hydrothermal.description)
- `hydrothermal_plastic_strain_reference::`[`ParameterFloat`](@ref): $(PDATA.hydrothermal_plastic_strain_reference.description)
- `sediment_thickness_threshold::`[`ParameterFloat`](@ref): $(PDATA.sediment_thickness_threshold.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
parameter objects

# Nested Dot Access
- `iuse_hydrothermal = $(ROOT_NAME).$(GRP_NAME).iuse_hydrothermal.value`
- `iuse_melt_lens = $(ROOT_NAME).$(GRP_NAME).iuse_melt_lens.value`
- `hydrothermal_smoothing_factor = $(ROOT_NAME).$(GRP_NAME).hydrothermal_smoothing_factor.value`
- `hydrothermal_nusselt_number = $(ROOT_NAME).$(GRP_NAME).hydrothermal_nusselt_number.value`
- `hydrothermal_max_temperature = $(ROOT_NAME).$(GRP_NAME).hydrothermal_max_temperature.value`
- `hydrothermal_max_depth = $(ROOT_NAME).$(GRP_NAME).hydrothermal_max_depth.value`
- `iuse_plastic_strain_rate_for_hydrothermal = $(ROOT_NAME).$(GRP_NAME).iuse_plastic_strain_rate_for_hydrothermal.value`
- `hydrothermal_decay_length = $(ROOT_NAME).$(GRP_NAME).hydrothermal_decay_length.value`
- `hydrothermal_buffer_distance = $(ROOT_NAME).$(GRP_NAME).hydrothermal_buffer_distance.value`
- `hydrothermal_plastic_strain_rate_reference = $(ROOT_NAME).$(GRP_NAME).hydrothermal_plastic_strain_rate_reference.value`
- `iuse_plastic_strain_for_hydrothermal = $(ROOT_NAME).$(GRP_NAME).iuse_plastic_strain_for_hydrothermal.value`
- `hydrothermal_plastic_strain_reference = $(ROOT_NAME).$(GRP_NAME).hydrothermal_plastic_strain_reference.value`
- `sediment_thickness_threshold = $(ROOT_NAME).$(GRP_NAME).sediment_thickness_threshold.value`

# Constructor
    Hydrothermal()

Create a new Hydrothermal parameter group with default values.

# Returns
- `Hydrothermal`: New Hydrothermal parameter group with initialized values

"""
mutable struct Hydrothermal <: AbstractParameterGroup
    iuse_hydrothermal::ParameterInt
    iuse_melt_lens::ParameterInt
    hydrothermal_smoothing_factor::ParameterFloat
    hydrothermal_nusselt_number::ParameterFloat
    hydrothermal_max_temperature::ParameterFloat
    hydrothermal_max_depth::ParameterFloat
    iuse_plastic_strain_rate_for_hydrothermal::ParameterInt
    hydrothermal_decay_length::ParameterFloat
    hydrothermal_buffer_distance::ParameterFloat
    hydrothermal_plastic_strain_rate_reference::ParameterFloat
    iuse_plastic_strain_for_hydrothermal::ParameterInt
    hydrothermal_plastic_strain_reference::ParameterFloat
    sediment_thickness_threshold::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Hydrothermal()::Hydrothermal
    pdata = get_eb_parameters()
    data = Hydrothermal(
        pdata.iuse_hydrothermal,
        pdata.iuse_melt_lens,
        pdata.hydrothermal_smoothing_factor,
        pdata.hydrothermal_nusselt_number,
        pdata.hydrothermal_max_temperature,
        pdata.hydrothermal_max_depth,
        pdata.iuse_plastic_strain_rate_for_hydrothermal,
        pdata.hydrothermal_decay_length,
        pdata.hydrothermal_buffer_distance,
        pdata.hydrothermal_plastic_strain_rate_reference,
        pdata.iuse_plastic_strain_for_hydrothermal,
        pdata.hydrothermal_plastic_strain_reference,
        pdata.sediment_thickness_threshold,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
