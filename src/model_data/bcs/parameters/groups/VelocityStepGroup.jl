"""
Module for velocity step parameters.

Provides data structures for configuring time-dependent velocity changes, 
allowing velocity to increase or decrease by specified factors at specified 
times.
"""
module VelocityStepGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: AbstractParameterGroup
import EarthBox.Parameters: ParameterFloat, ParameterInt

const ROOT_NAME = "model.bcs.parameters"
const GRP_NAME = "velocity_step"

const PDATA = get_eb_parameters()

"""
    VelocityStep <: AbstractParameterGroup

Parameter group for velocity step changes.

# Fields
- `iuse_velocity_step::`[`ParameterInt`](@ref): $(PDATA.iuse_velocity_step.description)
- `velocity_step_factor::`[`ParameterFloat`](@ref): $(PDATA.velocity_step_factor.description)
- `timestep_adjustment_factor::`[`ParameterFloat`](@ref): $(PDATA.timestep_adjustment_factor.description)
- `velocity_step_time::`[`ParameterFloat`](@ref): $(PDATA.velocity_step_time.description)
- `velocity_second_step_time::`[`ParameterFloat`](@ref): $(PDATA.velocity_second_step_time.description)
- `ivelocity_step_counter::`[`ParameterInt`](@ref): $(PDATA.ivelocity_step_counter.description)
- `velocity_second_step_factor::`[`ParameterFloat`](@ref): $(PDATA.velocity_second_step_factor.description)
- `timestep_second_adjustment_factor::`[`ParameterFloat`](@ref): $(PDATA.timestep_second_adjustment_factor.description)

# Nested Dot Access
- `iuse_velocity_step = $(ROOT_NAME).$(GRP_NAME).iuse_velocity_step.value`
- `velocity_step_factor = $(ROOT_NAME).$(GRP_NAME).velocity_step_factor.value`
- `timestep_adjustment_factor = $(ROOT_NAME).$(GRP_NAME).timestep_adjustment_factor.value`
- `velocity_step_time = $(ROOT_NAME).$(GRP_NAME).velocity_step_time.value`
- `velocity_second_step_time = $(ROOT_NAME).$(GRP_NAME).velocity_second_step_time.value`
- `ivelocity_step_counter = $(ROOT_NAME).$(GRP_NAME).ivelocity_step_counter.value`
- `velocity_second_step_factor = $(ROOT_NAME).$(GRP_NAME).velocity_second_step_factor.value`
- `timestep_second_adjustment_factor = $(ROOT_NAME).$(GRP_NAME).timestep_second_adjustment_factor.value`

# Constructor
    VelocityStep()

Initializes velocity step parameters with default values.
"""
mutable struct VelocityStep <: AbstractParameterGroup
    iuse_velocity_step::ParameterInt
    velocity_step_factor::ParameterFloat
    timestep_adjustment_factor::ParameterFloat
    velocity_step_time::ParameterFloat
    velocity_second_step_time::ParameterFloat
    ivelocity_step_counter::ParameterInt
    velocity_second_step_factor::ParameterFloat
    timestep_second_adjustment_factor::ParameterFloat
end

function VelocityStep()::VelocityStep
    pdata = get_eb_parameters()
    data = VelocityStep(
        pdata.iuse_velocity_step,
        pdata.velocity_step_factor,
        pdata.timestep_adjustment_factor,
        pdata.velocity_step_time,
        pdata.velocity_second_step_time,
        pdata.ivelocity_step_counter,
        pdata.velocity_second_step_factor,
        pdata.timestep_second_adjustment_factor
    )
    return data
end

end # module
