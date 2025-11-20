module CarbArray2D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray2D

"""
    CarbArray2DState

A struct for storing 2D carbonate data arrays with metadata.

# Fields
- `array::Matrix{Float64}`: The 2D array of float64 values
- `name::String`: Name of the array
- `units::Vector{String}`: Units for the array values
- `description::Vector{String}`: Description of the array contents
"""
mutable struct CarbArray2DState <: AbstractEarthBoxArray2D
    array::Matrix{Float64}
    name::String
    units::Vector{String}
    description::Vector{String}
end

end # module 