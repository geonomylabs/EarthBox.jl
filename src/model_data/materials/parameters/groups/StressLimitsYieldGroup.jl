module StressLimitsYieldGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "stress_limits_yield"

const PDATA = get_eb_parameters()

"""
    StressLimitsYield <: AbstractParameterGroup

Parameter group for yield stress limits and plastic failure properties.

# Fields
- `yield_stress_min::`[`ParameterFloat`](@ref): $(PDATA.yield_stress_min.description)
- `yield_stress_max::`[`ParameterFloat`](@ref): $(PDATA.yield_stress_max.description)
- `iuse_fluid_pressure_for_yield::`[`ParameterInt`](@ref): $(PDATA.iuse_fluid_pressure_for_yield.description)
- `plastic_healing_rate::`[`ParameterFloat`](@ref): $(PDATA.plastic_healing_rate.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `yield_stress_min = $(ROOT_NAME).$(GRP_NAME).yield_stress_min.value`
- `yield_stress_max = $(ROOT_NAME).$(GRP_NAME).yield_stress_max.value`
- `iuse_fluid_pressure_for_yield = $(ROOT_NAME).$(GRP_NAME).iuse_fluid_pressure_for_yield.value`
- `plastic_healing_rate = $(ROOT_NAME).$(GRP_NAME).plastic_healing_rate.value`

# Constructor
    StressLimitsYield()

Create a new StressLimitsYield parameter group with default values.

# Returns
- `StressLimitsYield`: New StressLimitsYield parameter group with initialized values

"""
mutable struct StressLimitsYield <: AbstractParameterGroup
    yield_stress_min::ParameterFloat
    yield_stress_max::ParameterFloat
    iuse_fluid_pressure_for_yield::ParameterInt
    plastic_healing_rate::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function StressLimitsYield()::StressLimitsYield
    pdata = get_eb_parameters()
    data = StressLimitsYield(
        pdata.yield_stress_min,
        pdata.yield_stress_max,
        pdata.iuse_fluid_pressure_for_yield,
        pdata.plastic_healing_rate,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
