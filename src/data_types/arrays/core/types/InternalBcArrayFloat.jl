"""
    InternalBcArrayFloat

Module containing the internal boundary condition float array type used throughout the EarthBox model.
This type represents a 1D array with associated metadata including name and descriptions
for internal boundary conditions.
"""
module InternalBcArrayFloat

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D

"""
    InternalBcArrayFloatState

Struct representing a 1D internal boundary condition float array with associated metadata.

# Fields
- `array`: the actual array data (Vector{Float64})
- `name`: identifier for the array
- `description`: array of detailed descriptions for each element
"""
mutable struct InternalBcArrayFloatState <: AbstractEarthBoxArray1D
    array::Vector{Float64}
    name::String
    description::Vector{String}
end

end # module 