module Grids2dContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ArrayCollection: Arrays
import .ParameterCollection: Parameters

export Grids

"""
    Grids <: CollectionContainer

Data structure containing parameter and array objects associated with 2D staggered grids.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for basic and staggered grids
- `arrays::`[`Arrays`](@ref): Array groups containing coordinates and spacing for basic and staggered 

# Constructor
    Grids(ynum::Int, xnum::Int, ysize::Float64, xsize::Float64)::Grids

Create a new Grids collection with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction
- `ysize::Float64`: Domain size in y-direction in meters
- `xsize::Float64`: Domain size in x-direction in meters

"""
mutable struct Grids <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function Grids(ynum::Int, xnum::Int, ysize::Float64, xsize::Float64)::Grids
    @assert ynum > 2 "Number of grid points in y-direction must be greater than 2"
    @assert xnum > 2 "Number of grid points in x-direction must be greater than 2"
    @assert ysize > 0.0 "Domain size in y-direction must be positive"
    @assert xsize > 0.0 "Domain size in x-direction must be positive"
    @assert ysize/ynum > 0.0 "Grid spacing in y-direction must be positive"
    @assert xsize/xnum > 0.0 "Grid spacing in x-direction must be positive"
    return Grids(Parameters(ynum, xnum, ysize, xsize), Arrays(ynum, xnum))
end

end # module 