"""
    InternalBcArrayInt

Module containing the internal boundary condition integer array type used throughout the EarthBox model.
This type represents a 1D array with associated metadata including name and descriptions
for internal boundary conditions.
"""
module InternalBcArrayInt

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D

"""
    InternalBcArrayIntState

Struct representing a 1D internal boundary condition integer array with associated metadata.

# Fields
- `array`: the actual array data (Vector{Int})
- `name`: identifier for the array
- `description`: array of detailed descriptions for each element
"""
mutable struct InternalBcArrayIntState <: AbstractEarthBoxArray1D
    array::Vector{Int}
    name::String
    description::Vector{String}
end

end # module 