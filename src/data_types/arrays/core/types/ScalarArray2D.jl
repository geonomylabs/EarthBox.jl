"""
    ScalarArray2D

Module for handling 2D scalar arrays in EarthBox. Provides functionality for creating and managing
2D arrays with different grid types (basic, pressure, vx, vy) and associated metadata.
"""
module ScalarArray2D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray2D
import ...ArrayUtils: OutputFormat

function grid_array2D(ynum::Int64, xnum::Int64, ::Val{:basic})::Array{Float64, 2}
    return zeros(Float64, ynum, xnum)
end

function grid_array2D(ynum::Int64, xnum::Int64, ::Val{:vx})::Array{Float64, 2}
    return zeros(Float64, ynum+1, xnum)
end

function grid_array2D(ynum::Int64, xnum::Int64, ::Val{:vy})::Array{Float64, 2}
    return zeros(Float64, ynum, xnum+1)
end

function grid_array2D(ynum::Int64, xnum::Int64, ::Val{:pressure})::Array{Float64, 2}
    return zeros(Float64, ynum-1, xnum-1)
end

"""
    ScalarArray2DState

Mutable struct representing a 2D scalar array with associated metadata.

# Fields
- `ny::Int64`: Number of nodes in the y-direction
- `nx::Int64`: Number of nodes in the x-direction
- `array::Matrix{Float64}`: The 2D array data
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `grid_type::String`: Type of grid ("basic", "pressure", "vx", "vy")
- `description::String`: Description of the array's purpose
- `outform::OutputFormat`: Output formatting specifications
"""
mutable struct ScalarArray2DState <: AbstractEarthBoxArray2D
    ny::Int64
    nx::Int64
    array::Matrix{Float64}
    name::String
    units::String
    grid_type::String
    description::String
    outform::OutputFormat
end

"""
    ScalarArray2DState(ny, nx, name, units, grid_type, description)

Construct a new ScalarArray2DState with the specified parameters.

# Arguments
- `ny::Int64`: Number of nodes in the y-direction
- `nx::Int64`: Number of nodes in the x-direction
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `grid_type::String`: Type of grid ("basic", "pressure", "vx", "vy")
- `description::String`: Description of the array's purpose

# Returns
- `ScalarArray2DState`: A new instance with initialized array and metadata

# Grid Types
- "basic": Standard grid with dimensions (ny, nx)
- "pressure": Pressure grid with dimensions (ny-1, nx-1)
- "vx": Velocity x-component grid with dimensions (ny+1, nx)
- "vy": Velocity y-component grid with dimensions (ny, nx+1)
"""
function ScalarArray2DState(
    ny::Int64,
    nx::Int64,
    name::String,
    units::String,
    grid_type::String,
    description::String
)::ScalarArray2DState
    if grid_type == "basic"
        array = grid_array2D(ny, nx, Val(:basic))
    elseif grid_type == "pressure"
        array = grid_array2D(ny, nx, Val(:pressure))
    elseif grid_type == "vx"
        array = grid_array2D(ny, nx, Val(:vx))
    elseif grid_type == "vy"
        array = grid_array2D(ny, nx, Val(:vy))
    end
    outform = OutputFormat(1.0, 0.0, units, name, false, "undefined.dat", name)
    return ScalarArray2DState(
        ny, nx, array, name, units, grid_type, description, outform
        )
end

end # module 