module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import EarthBox.Parameters: ParameterInt, ParameterStr
import .ThermalCondGroup: ThermalCond
import .HeatOptionsGroup: HeatOptions
import .TempChangeLimitGroup: TempChangeLimit
import .BuildGroup: Build
import .InitialConditionGroup: InitialCondition
import .SteadyStateGroup: SteadyState
import .RhoCpGroup: RhoCp

"""
    Parameters <: AbstractParameterCollection

Parameter collection for heat equation solver.

# Fields
- `thermalcond::`[`ThermalCond`](@ref): Thermal conductivity model options
- `heat_options::`[`HeatOptions`](@ref): Heat solver on/off switches
- `temp_change_limit::`[`TempChangeLimit`](@ref): Temperature change limits
- `build::`[`Build`](@ref): Sparse matrix build parameters
- `initial_condition::`[`InitialCondition`](@ref): Initial temperature conditions
- `steady_state::`[`SteadyState`](@ref): Steady state thermal structure parameters
- `rhocp::`[`RhoCp`](@ref): Density times heat capacity options

# Constructor
    Parameters(ynum::Int, xnum::Int)

Create a new Parameters collection with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

"""
mutable struct Parameters <: AbstractParameterCollection
    thermalcond::ThermalCond
    heat_options::HeatOptions
    temp_change_limit::TempChangeLimit
    build::Build
    initial_condition::InitialCondition
    steady_state::SteadyState
    rhocp::RhoCp
end

function Parameters(ynum::Int, xnum::Int)::Parameters
    return Parameters(
        ThermalCond(),
        HeatOptions(),
        TempChangeLimit(),
        Build(xnum, ynum),
        InitialCondition(),
        SteadyState(),
        RhoCp()
    )
end

end # module 