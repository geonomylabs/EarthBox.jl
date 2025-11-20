"""
    MarkerArrayFloat1D

Module for handling 1D float marker arrays in EarthBox. Provides functionality for 
creating and managing 1D float arrays associated with marker particles.
Supports various floating-point types including Float64, Float32, and Float16.
"""
module MarkerArrayFloat1D

import LinearAlgebra
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D
import ...ArrayUtils: OutputFormat

"""
    MarkerArrayFloat1DState{T<:AbstractFloat}

Mutable struct representing a 1D float marker array with associated metadata.

# Type Parameters
- `T<:AbstractFloat`: The floating-point type for the array (e.g., Float64, Float32, Float16)

# Fields
- `array::Vector{T}`: The 1D float array data for markers
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `description::String`: Description of the array's purpose
- `outform::OutputFormat`: Output formatting specifications
"""
mutable struct MarkerArrayFloat1DState{T<:AbstractFloat} <: AbstractEarthBoxArray1D
    array::Vector{T}
    name::String
    units::String
    description::String
    outform::OutputFormat
end

function MarkerArrayFloat1DState(
    marker_array::Vector{T},
    attr_name::String,
    units::String,
    description::String
)::MarkerArrayFloat1DState{T} where {T<:AbstractFloat}
    outform = OutputFormat(
        1.0,
        0.0,
        units,
        attr_name,
        false,
        "undefined.dat",
        attr_name
    )
    return MarkerArrayFloat1DState{T}(
        marker_array, attr_name, units, description, outform)
end

end # module
