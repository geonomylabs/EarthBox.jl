module TopographyContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters
import .ArrayCollection: Arrays

export Topography

"""
    Topography <: CollectionContainer

Data structure containing parameter and array objects for topography evolution.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for topography configuration
- `arrays::`[`Arrays`](@ref): Array groups for topographic grid data

# Constructor
    Topography()::Topography

Create a new Topography collection with default parameters.

"""
mutable struct Topography <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function Topography()::Topography
    parameters = Parameters()
    toponum = parameters.topo_grid.toponum.value
    arrays = Arrays(toponum)
    return Topography(parameters, arrays)
end

end # module 