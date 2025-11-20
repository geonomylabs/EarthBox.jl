module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .MarkerAndGridInterpOptionsGroup: MarkerAndGridInterpOptions

"""
    Parameters <: AbstractParameterCollection

Parameter collection for marker-grid interpolation.

# Fields
- `interp_options::`[`MarkerAndGridInterpOptions`](@ref): Interpolation method options

# Constructor
    Parameters()

Create a new Parameters collection with default interpolation options.

"""
mutable struct Parameters <: AbstractParameterCollection
    interp_options::MarkerAndGridInterpOptions
end

function Parameters()::Parameters
    return Parameters(MarkerAndGridInterpOptions())
end

end # module 