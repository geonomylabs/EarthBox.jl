module ConversionContainer

include("parameters/ParameterCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ParameterCollection: Parameters

"""
    Conversion

Parameter collection for conversion parameters.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for conversion configuration

"""
mutable struct Conversion <: CollectionContainer
    parameters::Parameters
end

function Conversion()::Conversion
    parameters = Parameters()
    return Conversion(parameters)
end

end # module 