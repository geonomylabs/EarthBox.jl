module GravityContainer

include("parameters/ParameterCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters

"""
    Gravity

Struct containing all gravity-related components for the EarthBox model.

# Fields
- `parameters`: Gravity parameters including gravity components and control flags

"""
mutable struct Gravity <: CollectionContainer
    parameters::Parameters
end

function Gravity()::Gravity
    parameters = Parameters()
    return Gravity(parameters)
end

end # module 