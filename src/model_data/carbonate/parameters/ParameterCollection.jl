module ParameterCollection

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import EarthBox.Parameters: ParameterInt, ParameterFloat

const ROOT_NAME = "model.carbonate.parameters"

const PDATA = get_eb_parameters()

"""
    Parameters <: AbstractParameterCollection

Parameter collection for carbonate deposition configuration.

# Fields
- `iuse_carb::`[`ParameterInt`](@ref): $(PDATA.iuse_carb.description)
- `photic_thick_m::`[`ParameterFloat`](@ref): $(PDATA.photic_thick_m.description)
- `carb_growth_rad::`[`ParameterFloat`](@ref): $(PDATA.carb_growth_rad.description)
- `carb_base_rate::`[`ParameterFloat`](@ref): $(PDATA.carb_base_rate.description)
- `carb_time_myr::`[`ParameterFloat`](@ref): $(PDATA.carb_time_myr.description)
- `carb_jump_time_myr::`[`ParameterFloat`](@ref): $(PDATA.carb_jump_time_myr.description)
- `carb_base_rate_jump::`[`ParameterFloat`](@ref): $(PDATA.carb_base_rate_jump.description)

# Nested Dot Access
- `iuse_carb = $(ROOT_NAME).iuse_carb.value`
- `photic_thick_m = $(ROOT_NAME).photic_thick_m.value`
- `carb_growth_rad = $(ROOT_NAME).carb_growth_rad.value`
- `carb_base_rate = $(ROOT_NAME).carb_base_rate.value`
- `carb_time_myr = $(ROOT_NAME).carb_time_myr.value`
- `carb_jump_time_myr = $(ROOT_NAME).carb_jump_time_myr.value`
- `carb_base_rate_jump = $(ROOT_NAME).carb_base_rate_jump.value`

# Constructor
    Parameters()

# Returns
- `Parameters`: New Parameters collection with initialized values

"""
mutable struct Parameters <: AbstractParameterCollection
    iuse_carb::ParameterInt
    photic_thick_m::ParameterFloat
    carb_growth_rad::ParameterFloat
    carb_base_rate::ParameterFloat
    carb_time_myr::ParameterFloat
    carb_jump_time_myr::ParameterFloat
    carb_base_rate_jump::ParameterFloat
end

function Parameters()::Parameters
    pdata = get_eb_parameters()
    return Parameters(
        pdata.iuse_carb,
        pdata.photic_thick_m,
        pdata.carb_growth_rad,
        pdata.carb_base_rate,
        pdata.carb_time_myr,
        pdata.carb_jump_time_myr,
        pdata.carb_base_rate_jump
    )
end

end # module
