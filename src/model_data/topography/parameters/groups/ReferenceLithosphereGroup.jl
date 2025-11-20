module ReferenceLithosphereGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.topography.parameters"
const GRP_NAME = "reference_lithosphere"

const PDATA = get_eb_parameters()

"""
    ReferenceLithosphere <: AbstractParameterGroup

Parameter group for reference lithosphere parameters.

# Fields
- `thickness_upper_continental_crust_ref::`[`ParameterFloat`](@ref): $(PDATA.thickness_upper_continental_crust_ref.description)
- `thickness_lower_continental_crust_ref::`[`ParameterFloat`](@ref): $(PDATA.thickness_lower_continental_crust_ref.description)
- `thickness_lithosphere_ref::`[`ParameterFloat`](@ref): $(PDATA.thickness_lithosphere_ref.description)
- `gridy_spacing_ref::`[`ParameterFloat`](@ref): $(PDATA.gridy_spacing_ref.description)
- `temperature_top_ref::`[`ParameterFloat`](@ref): $(PDATA.temperature_top_ref.description)
- `temperature_moho_ref::`[`ParameterFloat`](@ref): $(PDATA.temperature_moho_ref.description)
- `temperature_base_lith_ref::`[`ParameterFloat`](@ref): $(PDATA.temperature_base_lith_ref.description)
- `adiabatic_gradient_ref::`[`ParameterFloat`](@ref): $(PDATA.adiabatic_gradient_ref.description)
- `iuse_linear_segments::`[`ParameterInt`](@ref): $(PDATA.iuse_linear_segments.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `thickness_upper_continental_crust_ref = $(ROOT_NAME).$(GRP_NAME).thickness_upper_continental_crust_ref.value`
- `thickness_lower_continental_crust_ref = $(ROOT_NAME).$(GRP_NAME).thickness_lower_continental_crust_ref.value`
- `thickness_lithosphere_ref = $(ROOT_NAME).$(GRP_NAME).thickness_lithosphere_ref.value`
- `gridy_spacing_ref = $(ROOT_NAME).$(GRP_NAME).gridy_spacing_ref.value`
- `temperature_top_ref = $(ROOT_NAME).$(GRP_NAME).temperature_top_ref.value`
- `temperature_moho_ref = $(ROOT_NAME).$(GRP_NAME).temperature_moho_ref.value`
- `temperature_base_lith_ref = $(ROOT_NAME).$(GRP_NAME).temperature_base_lith_ref.value`
- `adiabatic_gradient_ref = $(ROOT_NAME).$(GRP_NAME).adiabatic_gradient_ref.value`
- `iuse_linear_segments = $(ROOT_NAME).$(GRP_NAME).iuse_linear_segments.value`

# Constructor
    ReferenceLithosphere()

# Returns
- `ReferenceLithosphere`: New ReferenceLithosphere parameter group with initialized values

"""
mutable struct ReferenceLithosphere <: AbstractParameterGroup
    thickness_upper_continental_crust_ref::ParameterFloat
    thickness_lower_continental_crust_ref::ParameterFloat
    thickness_lithosphere_ref::ParameterFloat
    gridy_spacing_ref::ParameterFloat
    temperature_top_ref::ParameterFloat
    temperature_moho_ref::ParameterFloat
    temperature_base_lith_ref::ParameterFloat
    adiabatic_gradient_ref::ParameterFloat
    iuse_linear_segments::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function ReferenceLithosphere()::ReferenceLithosphere
    pdata = get_eb_parameters()
    data = ReferenceLithosphere(
        pdata.thickness_upper_continental_crust_ref,
        pdata.thickness_lower_continental_crust_ref,
        pdata.thickness_lithosphere_ref,
        pdata.gridy_spacing_ref,
        pdata.temperature_top_ref,
        pdata.temperature_moho_ref,
        pdata.temperature_base_lith_ref,
        pdata.adiabatic_gradient_ref,
        pdata.iuse_linear_segments,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
