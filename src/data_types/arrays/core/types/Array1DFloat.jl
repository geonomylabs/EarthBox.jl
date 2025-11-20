"""
    Array1DFloat

Module containing a basic 1D float array type for the EarthBox model.
"""
module Array1DFloat

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D

"""
    Array1DFloatState

Mutable struct representing a 1D float array with associated metadata.

# Fields
- `array::Vector{Float64}`: The 1D float array data
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `description::String`: Description of the array's purpose
"""
mutable struct Array1DFloatState <: AbstractEarthBoxArray1D
    array::Vector{Float64}
    name::String
    units::String
    description::String
end

end # module