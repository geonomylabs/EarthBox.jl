"""
    BcArrayFloat

Module containing the boundary condition array type used throughout the EarthBox model.
This type represents a 2D array with associated metadata including name, units,
and description for boundary conditions.
"""
module BcArrayFloat

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray2D

"""
    BcArrayFloatState

Struct representing a 2D boundary condition array with associated metadata.

# Fields
- `array`: the actual array data (Matrix{Float64})
- `name`: identifier for the array
- `units`: physical units of the array values
- `description`: detailed description of what the array represents
"""
mutable struct BcArrayFloatState <: AbstractEarthBoxArray2D
    array::Matrix{Float64}
    name::String
    units::String
    description::String
end

end # module 