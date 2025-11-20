"""
    MaterialArrayInt2D

Module for handling 2D material arrays in EarthBox. Provides functionality for creating and managing
2D integer arrays for material properties with associated metadata.
"""
module MaterialArrayInt2D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray2D

"""
    MaterialArrayInt2DState

Mutable struct representing a 2D material array with associated metadata.

# Fields
- `array::Matrix{Int64}`: The 2D array data
- `name::String`: Name of the array
- `units::Vector{String}`: List of physical units
- `description::Vector{String}`: List of descriptions
"""
mutable struct MaterialArrayInt2DState <: AbstractEarthBoxArray2D
    array::Matrix{Int64}
    name::String
    units::Vector{String}
    description::Vector{String}
end

end # module 