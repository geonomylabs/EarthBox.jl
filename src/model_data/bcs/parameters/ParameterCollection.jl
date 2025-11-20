module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .BcOptionsGroup: BcOptions
import .TemperatureGroup: Temperature
import .VelocityGroup: Velocity
import .PressureGroup: Pressure
import .VelocityStepGroup: VelocityStep

"""
    Parameters

Collection of boundary condition parameters.

# Fields
- `bc_options::`[`BcOptions`](@ref): General boundary condition options
- `temperature::`[`Temperature`](@ref): Temperature boundary conditions
- `velocity::`[`Velocity`](@ref): Velocity boundary conditions  
- `pressure::`[`Pressure`](@ref): Pressure boundary conditions
- `velocity_step::`[`VelocityStep`](@ref): Velocity step boundary conditions
"""
mutable struct Parameters <: AbstractParameterCollection
    bc_options::BcOptions
    temperature::Temperature
    velocity::Velocity
    pressure::Pressure
    velocity_step::VelocityStep
end

function Parameters()::Parameters
    return Parameters(
        BcOptions(),
        Temperature(),
        Velocity(),
        Pressure(),
        VelocityStep()
    )
end

end # module 