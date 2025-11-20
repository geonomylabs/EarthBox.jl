"""
    GridArray1D

Module containing the basic 1D grid array type used throughout the EarthBox model.
This type represents a single array with associated metadata including name, units,
and description.
"""
module GridArray1D

import LinearAlgebra
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D

"""
    GridArray1DState

Struct representing a 1D grid array with associated metadata.

# Fields
- `array`: the actual array data (Vector{Float64})
- `name`: identifier for the array
- `units`: physical units of the array values
- `description`: detailed description of what the array represents
"""
mutable struct GridArray1DState <: AbstractEarthBoxArray1D
    array::Vector{Float64}
    name::String
    units::String
    description::String
end

end # module