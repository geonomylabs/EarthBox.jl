module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .GridWeightsGroup: GridWeights
import .MarkerWeightsGroup: MarkerWeights

"""
    Arrays <: AbstractArrayCollection

Data structure containing array groups for marker-grid interpolation weights.

# Fields
- `grid_weights::`[`GridWeights`](@ref): Summed weight arrays on grid nodes
- `marker_weights::`[`MarkerWeights`](@ref): Individual marker weight arrays

# Constructor
    Arrays(ynum::Int, xnum::Int, marknum::Int)::Arrays

Create a new Arrays collection with the given grid and marker dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction
- `marknum::Int`: Number of markers in the model

"""
mutable struct Arrays <: AbstractArrayCollection
    grid_weights::GridWeights
    marker_weights::MarkerWeights
end

function Arrays(ynum::Int, xnum::Int, marknum::Int)::Arrays
    return Arrays(GridWeights(ynum, xnum), MarkerWeights(marknum))
end

end # module 