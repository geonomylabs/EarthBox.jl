module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .MainTimeLoopGroup: MainTimeLoop
import .ThermalLoopGroup: ThermalLoop
import .OutputStepsGroup: OutputSteps
import .SingleIncreaseGroup: SingleIncrease
import .MultipleIncreaseGroup: MultipleIncrease
import .BoundaryDisplacementStoppingGroup: BoundaryDisplacementStopping

"""
    Parameters <: AbstractParameterCollection

Collection of timestep parameters.

# Fields
- `main_time_loop::`[`MainTimeLoop`](@ref): Parameters used for main model time loop
- `thermal_loop::`[`ThermalLoop`](@ref): Parameters used for adaptive time stepping loop for heat 
    equation
- `output_steps::`[`OutputSteps`](@ref): Parameters used for output steps
- `single_increase::`[`SingleIncrease`](@ref): Parameters used for modifying the time step by a 
    single step
- `multiple_increase::`[`MultipleIncrease`](@ref): Parameters used for modifying the time step over 
    multiple steps
- `boundary_displacement_stopping::`[`BoundaryDisplacementStopping`](@ref): Parameters used for 
    stopping the model when the displacement or extensional strain reaches a specified value
"""
mutable struct Parameters <: AbstractParameterCollection
    main_time_loop::MainTimeLoop
    thermal_loop::ThermalLoop
    output_steps::OutputSteps
    single_increase::SingleIncrease
    multiple_increase::MultipleIncrease
    boundary_displacement_stopping::BoundaryDisplacementStopping
end

function Parameters()::Parameters
    return Parameters(
        MainTimeLoop(),
        ThermalLoop(),
        OutputSteps(),
        SingleIncrease(),
        MultipleIncrease(),
        BoundaryDisplacementStopping()
    )
end

end # module 