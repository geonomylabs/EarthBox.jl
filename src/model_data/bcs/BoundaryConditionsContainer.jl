module BoundaryConditionsContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters
import .ArrayCollection: Arrays

export BoundaryConditions

"""
    BoundaryConditionsDataState

Struct containing boundary conditions parameters and arrays.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for boundary conditions
- `arrays::`[`Arrays`](@ref): Array groups for boundary conditions

"""
mutable struct BoundaryConditions <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function BoundaryConditions(ynum::Int, xnum::Int)::BoundaryConditions
    parameters = Parameters()
    arrays = Arrays(ynum, xnum)
    return BoundaryConditions(parameters, arrays)
end

end # module 