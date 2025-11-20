"""
    MaterialArrayInt1D

Module for handling 1D material arrays in EarthBox. Provides functionality for creating and managing
1D integer arrays for material properties with associated metadata.
"""
module MaterialArrayInt1D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D

"""
    MaterialArrayInt1DState

Mutable struct representing a 1D material array with associated metadata.

# Fields
- `array::Vector{Int64}`: The 1D array data
- `name::String`: Name of the array
- `units::Vector{String}`: List of physical units
- `description::Vector{String}`: List of descriptions
"""
mutable struct MaterialArrayInt1DState <: AbstractEarthBoxArray1D
    array::Vector{Int64}
    name::String
    units::Vector{String}
    description::Vector{String}
end

end # module 