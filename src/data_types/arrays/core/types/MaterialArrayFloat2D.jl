"""
    MaterialArrayFloat2D

Module for handling 2D material arrays in EarthBox. Provides functionality for creating and managing
2D floating-point arrays for material properties with associated metadata.
"""
module MaterialArrayFloat2D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray2D

"""
    MaterialArrayFloat2DState

Mutable struct representing a 2D material array with associated metadata.

# Fields
- `array::Matrix{Float64}`: The 2D array data
- `name::String`: Name of the array
- `units::Vector{String}`: List of physical units
- `description::Vector{String}`: List of descriptions
"""
mutable struct MaterialArrayFloat2DState <: AbstractEarthBoxArray2D
    array::Matrix{Float64}
    name::String
    units::Vector{String}
    description::Vector{String}
end

end # module 