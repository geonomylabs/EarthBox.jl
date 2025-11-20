"""
    SolutionArray1D

Module for handling 1D solution arrays in EarthBox. Provides functionality for creating and managing
1D arrays that represent solutions to various problems, with support for different solution types
(normal and velocity).
"""
module SolutionArray1D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D
import ...ArrayUtils: OutputFormat

"""
    SolutionArray1DState

Mutable struct representing a 1D solution array with associated metadata.

# Fields
- `array::Vector{Float64}`: The 1D array data
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `grid_type::String`: Type of grid (always "solution")
- `description::String`: Description of the array's purpose
- `outform::OutputFormat`: Output formatting specifications
"""
mutable struct SolutionArray1DState <: AbstractEarthBoxArray1D
    array::Vector{Float64}
    name::String
    units::String
    grid_type::String
    description::String
    outform::OutputFormat
end

"""
    SolutionArray1DState(ynum, xnum, name, units, description, solution_type)

Construct a new SolutionArray1DState with the specified parameters.

# Arguments
- `ynum::Int64`: Number of nodes in the y-direction
- `xnum::Int64`: Number of nodes in the x-direction
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `description::String`: Description of the array's purpose
- `solution_type::String`: Type of solution ("normal" or "velocity")

# Returns
- `SolutionArray1DState`: A new instance with initialized array and metadata

# Solution Types
- "normal": Creates an array of size (xnum - 1) * (ynum - 1) * 3
- "velocity": Creates an array of size (ynum + 1) * xnum + (xnum + 1) * ynum
"""
function SolutionArray1DState(
    ynum::Int64,
    xnum::Int64,
    name::String,
    units::String,
    description::String,
    solution_type::String
)::SolutionArray1DState
    if solution_type == "normal"
        nall = (xnum - 1) * (ynum - 1) * 3
    elseif solution_type == "velocity"
        nall = (ynum + 1) * xnum + (xnum + 1) * ynum
    end
    return SolutionArray1DState(
        zeros(Float64, nall),
        name,
        units,
        "solution",
        description,
        OutputFormat(1.0, 0.0, units, name, false, "undefined.dat", name)
    )
end

end # module 