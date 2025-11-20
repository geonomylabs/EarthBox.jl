"""
    RhsHeatArray1D

Module containing 1D array structure for right-hand side values in heat equations.
"""
module RhsHeatArray1D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D

"""
    RhsHeatArray1DState

Struct containing 1D array for heat equation right-hand side.

# Fields
- `array`: 1D array of float64 values
- `name`: Name of the array
- `units`: Units of the array values
- `description`: Description of what the array represents
"""
mutable struct RhsHeatArray1DState <: AbstractEarthBoxArray1D
    array::Vector{Float64}
    name::String
    units::String
    description::String
end

"""
    RhsHeatArray1DState(ynum::Int, xnum::Int, name::String, units::String, description::String)::RhsHeatArray1DState

Construct a new RhsHeatArray1DState with arrays of appropriate size.

# Arguments
- `ynum`: number of nodes in y-direction
- `xnum`: number of nodes in x-direction
- `name`: name of the array
- `units`: units of the array values
- `description`: description of what the array represents

# Returns
A new RhsHeatArray1DState with initialized arrays
"""
function RhsHeatArray1DState(
    ynum::Int, 
    xnum::Int, 
    name::String, 
    units::String, 
    description::String
)::RhsHeatArray1DState
    array = zeros(Float64, xnum*ynum)
    return RhsHeatArray1DState(array, name, units, description)
end

end # module 