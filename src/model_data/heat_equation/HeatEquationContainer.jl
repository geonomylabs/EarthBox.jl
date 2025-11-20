module HeatEquationContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import ..Grids2dContainer: Grids
import .ParameterCollection: Parameters
import .ArrayCollection: Arrays

export HeatEquation

"""
    HeatEquation <: CollectionContainer

Data structure containing parameter and array objects for heat equation solver.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for heat solver configuration
- `arrays::`[`Arrays`](@ref): Array groups for temperature and thermal properties

# Constructor
    HeatEquation(ynum::Int, xnum::Int)::HeatEquation

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

"""
mutable struct HeatEquation <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function HeatEquation(ynum::Int, xnum::Int)::HeatEquation
    @assert ynum > 0 "Number of grid points in y-direction must be positive"
    @assert xnum > 0 "Number of grid points in x-direction must be positive"
    return HeatEquation(Parameters(ynum, xnum), Arrays(ynum, xnum))
end

end # module 