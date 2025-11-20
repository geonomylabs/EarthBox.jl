"""
Module for temperature boundary condition parameters.

Provides data structures for configuring temperature boundary conditions along
domain boundaries and transient temperature conditions.
"""
module TemperatureGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: AbstractParameterGroup
import EarthBox.Parameters: ParameterFloat, ParameterInt

const ROOT_NAME = "model.bcs.parameters"
const GRP_NAME = "temperature"

const PDATA = get_eb_parameters()

"""
    Temperature <: AbstractParameterGroup

Parameter group for temperature boundary conditions.

# Fields
- `temperature_top::`[`ParameterFloat`](@ref): $(PDATA.temperature_top.description)
- `temperature_bottom::`[`ParameterFloat`](@ref): $(PDATA.temperature_bottom.description)
- `temperature_left::`[`ParameterFloat`](@ref): $(PDATA.temperature_left.description)
- `temperature_right::`[`ParameterFloat`](@ref): $(PDATA.temperature_right.description)
- `iuse_bottom_transient::`[`ParameterInt`](@ref): $(PDATA.iuse_bottom_transient.description)
- `temperature_bottom_transient::`[`ParameterFloat`](@ref): $(PDATA.temperature_bottom_transient.description)
- `temperature_bottom_original::`[`ParameterFloat`](@ref): $(PDATA.temperature_bottom_original.description)
- `start_time_bottom_transient::`[`ParameterFloat`](@ref): $(PDATA.start_time_bottom_transient.description)
- `end_time_bottom_transient::`[`ParameterFloat`](@ref): $(PDATA.end_time_bottom_transient.description)
- `delta_temperature_transient::`[`ParameterFloat`](@ref): $(PDATA.delta_temperature_transient.description)
- `temperature_base_lith_warmer_initial::`[`ParameterFloat`](@ref): $(PDATA.temperature_base_lith_warmer_initial.description)
- `temperature_bottom_warmer_initial::`[`ParameterFloat`](@ref): $(PDATA.temperature_bottom_warmer_initial.description)
- `temperature_bottom_cooler_final::`[`ParameterFloat`](@ref): $(PDATA.temperature_bottom_cooler_final.description)

# Nested Dot Access
- `temperature_top = $(ROOT_NAME).$(GRP_NAME).temperature_top.value`
- `temperature_bottom = $(ROOT_NAME).$(GRP_NAME).temperature_bottom.value`
- `temperature_left = $(ROOT_NAME).$(GRP_NAME).temperature_left.value`
- `temperature_right = $(ROOT_NAME).$(GRP_NAME).temperature_right.value`
- `iuse_bottom_transient = $(ROOT_NAME).$(GRP_NAME).iuse_bottom_transient.value`
- `temperature_bottom_transient = $(ROOT_NAME).$(GRP_NAME).temperature_bottom_transient.value`
- `temperature_bottom_original = $(ROOT_NAME).$(GRP_NAME).temperature_bottom_original.value`
- `start_time_bottom_transient = $(ROOT_NAME).$(GRP_NAME).start_time_bottom_transient.value`
- `end_time_bottom_transient = $(ROOT_NAME).$(GRP_NAME).end_time_bottom_transient.value`
- `delta_temperature_transient = $(ROOT_NAME).$(GRP_NAME).delta_temperature_transient.value`
- `temperature_base_lith_warmer_initial = $(ROOT_NAME).$(GRP_NAME).temperature_base_lith_warmer_initial.value`
- `temperature_bottom_warmer_initial = $(ROOT_NAME).$(GRP_NAME).temperature_bottom_warmer_initial.value`
- `temperature_bottom_cooler_final = $(ROOT_NAME).$(GRP_NAME).temperature_bottom_cooler_final.value`

# Constructor
    Temperature()

Initializes temperature boundary condition parameters with default values.
"""
mutable struct Temperature <: AbstractParameterGroup
    temperature_top::ParameterFloat
    temperature_bottom::ParameterFloat
    temperature_left::ParameterFloat
    temperature_right::ParameterFloat
    iuse_bottom_transient::ParameterInt
    temperature_bottom_transient::ParameterFloat
    temperature_bottom_original::ParameterFloat
    start_time_bottom_transient::ParameterFloat
    end_time_bottom_transient::ParameterFloat
    delta_temperature_transient::ParameterFloat
    temperature_base_lith_warmer_initial::ParameterFloat
    temperature_bottom_warmer_initial::ParameterFloat
    temperature_bottom_cooler_final::ParameterFloat
end

function Temperature()::Temperature
    pdata = get_eb_parameters()
    data = Temperature(
        pdata.temperature_top,
        pdata.temperature_bottom,
        pdata.temperature_left,
        pdata.temperature_right,
        pdata.iuse_bottom_transient,
        pdata.temperature_bottom_transient,
        pdata.temperature_bottom_original,
        pdata.start_time_bottom_transient,
        pdata.end_time_bottom_transient,
        pdata.delta_temperature_transient,
        pdata.temperature_base_lith_warmer_initial,
        pdata.temperature_bottom_warmer_initial,
        pdata.temperature_bottom_cooler_final
    )
    return data
end

end # module
