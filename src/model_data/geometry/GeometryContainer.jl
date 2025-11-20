module GeometryContainer

include("parameters/ParameterCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters
import ..DocTools: make_collection_doc

export Geometry

"""
    Geometry <: CollectionContainer

Struct containing material geometry parameters.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for material geometry.

"""
mutable struct Geometry <: CollectionContainer
    parameters::Parameters
end

function Geometry()::Geometry
    parameters = Parameters()
    return Geometry(parameters)
end

end # module 