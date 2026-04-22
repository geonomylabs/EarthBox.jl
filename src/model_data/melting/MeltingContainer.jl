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
    Melting(marknum::Int)::Melting

Create a new Melting collection with default values. `marknum` sizes the
marker-length scratch buffers in `arrays.buffers`.

"""
mutable struct Melting <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function Melting(marknum::Int)::Melting
    return Melting(Parameters(), Arrays(marknum))
end

end # module 