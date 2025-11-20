module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import EarthBox.Parameters: ParameterInt, ParameterStr
import .PicardGroup: Picard
import .SolutionNormsStokesGroup: SolutionNormsStokes
import .ResidualNormsStokesGroup: ResidualNormsStokes
import .VelocityCalcOptionsGroup: VelocityCalcOptions
import .OutputStokesGroup: OutputStokes
import .BuildSysStokesGroup: BuildSysStokes

"""
    Parameters <: AbstractParameterCollection

Parameter collection for Stokes-continuity solver.

# Fields
- `itype_density::`[`ParameterInt`](@ref): Density model ID
- `stype_density::`[`ParameterStr`](@ref): Density model name
- `picard::`[`Picard`](@ref): Picard iteration parameters
- `solution_norms::`[`SolutionNormsStokes`](@ref): Solution norm tracking
- `residual_norms::`[`ResidualNormsStokes`](@ref): Residual norm tracking
- `velocity_calc_option::`[`VelocityCalcOptions`](@ref): Velocity calculation options
- `output::`[`OutputStokes`](@ref): Output parameters
- `build::`[`BuildSysStokes`](@ref): Sparse matrix build parameters

# Constructor
    Parameters(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

"""
mutable struct Parameters <: AbstractParameterCollection
    itype_density::ParameterInt
    stype_density::ParameterStr
    picard::Picard
    solution_norms::SolutionNormsStokes
    residual_norms::ResidualNormsStokes
    velocity_calc_option::VelocityCalcOptions
    output::OutputStokes
    build::BuildSysStokes
end

function Parameters(ynum::Int, xnum::Int)::Parameters
    return Parameters(
        ParameterInt(
            0,                 # value
            "itype_density",   # name
            "None",            # units
            "Integer flag for density model."
        ),
        ParameterStr(
            "None",            # value
            "stype_density",   # name
            "None",            # units
            "Density model option name."
        ),
        Picard(),
        SolutionNormsStokes(),
        ResidualNormsStokes(),
        VelocityCalcOptions(),
        OutputStokes(),
        BuildSysStokes(ynum, xnum)
    )
end

end # module 