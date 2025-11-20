"""
    MaterialArrayFloat1D

Module for handling 1D material arrays in EarthBox. Provides functionality for creating and managing
1D floating-point arrays for material properties with associated metadata.
"""
module MaterialArrayFloat1D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D

"""
    MaterialArrayFloat1DState

Mutable struct representing a 1D material array with associated metadata.

# Fields
- `array::Vector{Float64}`: The 1D array data
- `name::String`: Name of the array
- `units::Vector{String}`: List of physical units
- `description::Vector{String}`: List of descriptions
"""
mutable struct MaterialArrayFloat1DState <: AbstractEarthBoxArray1D
    array::Vector{Float64}
    name::String
    units::Vector{String}
    description::Vector{String}
end

end # module 