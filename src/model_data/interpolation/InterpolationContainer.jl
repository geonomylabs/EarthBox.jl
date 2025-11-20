module InterpolationContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import ..Grids2dContainer: Grids
import ..MarkerContainer: Markers
import .ParameterCollection: Parameters
import .ArrayCollection: Arrays

export Interpolation

"""
    Interpolation <: CollectionContainer

Data structure containing parameter and array objects for marker-grid interpolation.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for interpolation options
- `arrays::`[`Arrays`](@ref): Array groups for grid and marker interpolation weights

# Constructor
    Interpolation(ynum::Int, xnum::Int, marknum::Int)::Interpolation

Create a new Interpolation collection with the given grid and marker dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction
- `marknum::Int`: Number of markers in the model

"""
mutable struct Interpolation <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function Interpolation(ynum::Int, xnum::Int, marknum::Int)::Interpolation
    @assert ynum > 0 "Number of grid points in y-direction must be positive"
    @assert xnum > 0 "Number of grid points in x-direction must be positive"
    @assert marknum > 0 "Number of markers must be positive"
    return Interpolation(Parameters(), Arrays(ynum, xnum, marknum))
end

function clear_grid_weights!(interp::Interpolation)
    interp.arrays.grid_weights.wtnodes.array .= 0.0
    interp.arrays.grid_weights.wtetas.array .= 0.0
    interp.arrays.grid_weights.wtetan.array .= 0.0
    interp.arrays.grid_weights.wtnodes_vy.array .= 0.0
end

end # module 