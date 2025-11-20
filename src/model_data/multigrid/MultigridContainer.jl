module MultiGrids2dContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ArrayCollection: Arrays
import .ParameterCollection: Parameters

const MAX_LEVELS = 10

mutable struct MultigridData <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end
    
"""
    MultigridData(nvcycles::Int64, levelnum::Int64, smoothing_iterations_on_finest_level::Int64)

# Arguments
- `nvcycles::Int64`: Number of multigrid iteration cycles.
- `levelnum::Int64`: Number of multigrid levels.
- `smoothing_iterations_on_finest_level::Int64`: Number of Gauss-Seidel smoothing iterations on the finest level.
"""
function MultigridData(;
    nvcycles::Int64, 
    levelnum::Int64, 
    smoothing_iterations_on_finest_level::Int64
)::MultigridData
    if levelnum > MAX_LEVELS
        error("Maximum number of multigrid levels is $MAX_LEVELS")
    end
    return MultigridData(
        Parameters(nvcycles, levelnum, MAX_LEVELS),
        Arrays(nvcycles, smoothing_iterations_on_finest_level, MAX_LEVELS)
    )
end

end