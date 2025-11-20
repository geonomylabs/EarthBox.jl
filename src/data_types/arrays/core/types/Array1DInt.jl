"""
    Array1DInt

Module containing a basic 1D integer array type for the EarthBox model.
"""
module Array1DInt

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D

"""
    Array1DIntState

Mutable struct representing a 1D integer array with associated metadata.

# Fields
- `array::Vector{Int}`: The 1D integer array data
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `description::String`: Description of the array's purpose
"""
mutable struct Array1DIntState <: AbstractEarthBoxArray1D
    array::Vector{Int}
    name::String
    units::String
    description::String
end

end # module