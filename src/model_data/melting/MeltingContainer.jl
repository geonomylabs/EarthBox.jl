module MeltingContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ArrayCollection: Arrays
import .ParameterCollection: Parameters

export Melting

"""
    Melting <: CollectionContainer

Data structure containing parameter and array objects for melting and melt extraction.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for melting options and configuration
- `arrays::`[`Arrays`](@ref): Array groups for melt extraction tracking

# Constructor
    Melting()::Melting

Create a new Melting collection with default values.

"""
mutable struct Melting <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function Melting()::Melting
    return Melting(Parameters(), Arrays())
end

end # module 