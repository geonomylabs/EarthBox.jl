module DepoAndErosionRatesGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.topography.parameters"
const GRP_NAME = "depo_and_erosion_rates"

const PDATA = get_eb_parameters()

"""
    DepoAndErosionRates <: AbstractParameterGroup

Parameter group for deposition and erosion rates.

# Fields
- `erosion_rate::`[`ParameterFloat`](@ref): $(PDATA.erosion_rate.description)
- `sedimentation_rate::`[`ParameterFloat`](@ref): $(PDATA.sedimentation_rate.description)
- `pelagic_sedimentation_rate::`[`ParameterFloat`](@ref): $(PDATA.pelagic_sedimentation_rate.description)
- `pelagic_sedimentation_rate_reduction_factor::`[`ParameterFloat`](@ref): $(PDATA.pelagic_sedimentation_rate_reduction_factor.description)
- `pelagic_sedimentation_rate_reduction_time::`[`ParameterFloat`](@ref): $(PDATA.pelagic_sedimentation_rate_reduction_time.description)
- `salt_deposition_rate::`[`ParameterFloat`](@ref): $(PDATA.salt_deposition_rate.description)
- `salt_start_time::`[`ParameterFloat`](@ref): $(PDATA.salt_start_time.description)
- `salt_end_time::`[`ParameterFloat`](@ref): $(PDATA.salt_end_time.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `erosion_rate = $(ROOT_NAME).$(GRP_NAME).erosion_rate.value`
- `sedimentation_rate = $(ROOT_NAME).$(GRP_NAME).sedimentation_rate.value`
- `pelagic_sedimentation_rate = $(ROOT_NAME).$(GRP_NAME).pelagic_sedimentation_rate.value`
- `pelagic_sedimentation_rate_reduction_factor = $(ROOT_NAME).$(GRP_NAME).pelagic_sedimentation_rate_reduction_factor.value`
- `pelagic_sedimentation_rate_reduction_time = $(ROOT_NAME).$(GRP_NAME).pelagic_sedimentation_rate_reduction_time.value`
- `salt_deposition_rate = $(ROOT_NAME).$(GRP_NAME).salt_deposition_rate.value`
- `salt_start_time = $(ROOT_NAME).$(GRP_NAME).salt_start_time.value`
- `salt_end_time = $(ROOT_NAME).$(GRP_NAME).salt_end_time.value`

# Constructor
    DepoAndErosionRates()

# Returns
- `DepoAndErosionRates`: New DepoAndErosionRates parameter group with initialized values

"""
mutable struct DepoAndErosionRates <: AbstractParameterGroup
    erosion_rate::ParameterFloat
    sedimentation_rate::ParameterFloat
    pelagic_sedimentation_rate::ParameterFloat
    pelagic_sedimentation_rate_reduction_factor::ParameterFloat
    pelagic_sedimentation_rate_reduction_time::ParameterFloat
    salt_deposition_rate::ParameterFloat
    salt_start_time::ParameterFloat
    salt_end_time::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function DepoAndErosionRates()::DepoAndErosionRates
    pdata = get_eb_parameters()
    data = DepoAndErosionRates(
        pdata.erosion_rate,
        pdata.sedimentation_rate,
        pdata.pelagic_sedimentation_rate,
        pdata.pelagic_sedimentation_rate_reduction_factor,
        pdata.pelagic_sedimentation_rate_reduction_time,
        pdata.salt_deposition_rate,
        pdata.salt_start_time,
        pdata.salt_end_time,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
