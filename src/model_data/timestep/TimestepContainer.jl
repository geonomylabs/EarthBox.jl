module TimestepContainer

include("parameters/ParameterCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters

export Timestep

"""
    Timestep <: CollectionContainer

Data structure containing parameter objects associated with model time stepping.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for model time stepping

"""
mutable struct Timestep <: CollectionContainer
    parameters::Parameters
end

function Timestep()::Timestep
    parameters = Parameters()
    return Timestep(parameters)
end

end # module 