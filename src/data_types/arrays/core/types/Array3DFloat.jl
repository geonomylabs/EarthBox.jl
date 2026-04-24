"""
    Array3DFloat

Module containing a basic 3D Float64 array type for the EarthBox model.
"""
module Array3DFloat

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray3D

"""
    Array3DFloatState

Mutable struct representing a 3D Float64 array with associated metadata.

# Fields
- `array::Array{Float64, 3}`: The 3D array data
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `description::String`: Description of the array's purpose
"""
mutable struct Array3DFloatState <: AbstractEarthBoxArray3D
    array::Array{Float64, 3}
    name::String
    units::String
    description::String
end

end # module
