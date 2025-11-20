module CarbonateContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters
import .ArrayCollection: Arrays

"""
    Carbonate <: CollectionContainer

Collection container for carbonate deposition configuration.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for carbonate deposition configuration
- `arrays::`[`Arrays`](@ref): Carbonate arrays containing grid nodes and properties
    
# Constructor
    Carbonate(toponum::Int64)::Carbonate

# Returns
- `Carbonate`: New Carbonate collection with initialized values

"""
mutable struct Carbonate <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function Carbonate(toponum::Int64)::Carbonate
    parameters = Parameters()
    arrays = Arrays(toponum)
    return Carbonate(parameters, arrays)
end

end # module 