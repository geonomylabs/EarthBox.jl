"""
    MarkerArrayInt1D

Module for handling 1D integer marker arrays in EarthBox. Provides functionality for 
creating and managing 1D integer arrays associated with marker particles.
Supports both Int64 and Int32 types.
"""
module MarkerArrayInt1D

import LinearAlgebra
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D
import ...ArrayUtils: OutputFormat

"""
    MarkerArrayInt1DState{T<:Integer}

Mutable struct representing a 1D integer marker array with associated metadata.

# Type Parameters
- `T<:Integer`: The integer type for the array (e.g., Int64, Int32, Int16)

# Fields
- `array::Vector{T}`: The 1D integer array data for markers
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `description::String`: Description of the array's purpose
- `outform::OutputFormat`: Output formatting specifications
"""
mutable struct MarkerArrayInt1DState{T<:Integer} <: AbstractEarthBoxArray1D
    array::Vector{T}
    name::String
    units::String
    description::String
    outform::OutputFormat
end

function MarkerArrayInt1DState(
    marker_array::Vector{T},
    attr_name::String,
    units::String,
    description::String
)::MarkerArrayInt1DState{T} where {T<:Integer}
    outform = OutputFormat(
        1.0,
        0.0,
        units,
        attr_name,
        false,
        "undefined.dat",
        attr_name
    )
    return MarkerArrayInt1DState{T}(
        marker_array, attr_name, units, description, outform)
end

end # end module
