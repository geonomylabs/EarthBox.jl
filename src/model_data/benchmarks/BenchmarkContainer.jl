module BenchmarkContainer

include("parameters/ParameterCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters

"""
    Benchmarks <: CollectionContainer

Parameter collection for benchmark parameters.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for benchmark configuration

# Constructor
    Benchmarks()

# Returns
- `Benchmarks`: New Benchmarks parameter collection with initialized values

"""
mutable struct Benchmarks <: CollectionContainer
    parameters::Parameters
end

function Benchmarks()::Benchmarks
    parameters = Parameters()
    return Benchmarks(parameters)
end

end # module 