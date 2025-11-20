module StokesContinuityContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters
import .ArrayCollection: Arrays

export StokesContinuity

"""
    StokesContinuity <: CollectionContainer

Data structure containing parameter and array objects for Stokes-continuity solver.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for Stokes solver configuration
- `arrays::`[`Arrays`](@ref): Array groups for velocity, pressure, and stress fields

# Constructor
    StokesContinuity(ynum::Int, xnum::Int)::StokesContinuity

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

"""
mutable struct StokesContinuity <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function StokesContinuity(ynum::Int, xnum::Int)::StokesContinuity
    @assert ynum > 0 "Number of grid points in y-direction must be positive"
    @assert xnum > 0 "Number of grid points in x-direction must be positive"
    return StokesContinuity(Parameters(ynum, xnum), Arrays(ynum, xnum))
end

end # module 