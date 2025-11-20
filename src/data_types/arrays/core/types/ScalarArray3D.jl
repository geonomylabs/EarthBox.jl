module ScalarArray3D

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray3D
import ...ArrayUtils: OutputFormat

function grid_array3D(ynum::Int64, xnum::Int64, znum::Int64, ::Val{:basic})::Array{Float64, 3}
    return zeros(Float64, ynum, xnum, znum)
end

function grid_array3D(ynum::Int64, xnum::Int64, znum::Int64, ::Val{:vx})::Array{Float64, 3}
    return zeros(Float64, ynum+1, xnum, znum+1)
end

function grid_array3D(ynum::Int64, xnum::Int64, znum::Int64, ::Val{:vy})::Array{Float64, 3}
    return zeros(Float64, ynum, xnum+1, znum+1)
end

function grid_array3D(ynum::Int64, xnum::Int64, znum::Int64, ::Val{:vz})::Array{Float64, 3}
    return zeros(Float64, ynum+1, xnum+1, znum)
end

function grid_array3D(ynum::Int64, xnum::Int64, znum::Int64, ::Val{:pressure})::Array{Float64, 3}
    return zeros(Float64, ynum-1, xnum-1, znum-1)
end

function grid_array3D(ynum::Int64, xnum::Int64, znum::Int64, ::Val{:shearxy})::Array{Float64, 3}
    return zeros(Float64, ynum, xnum, znum-1)
end

function grid_array3D(ynum::Int64, xnum::Int64, znum::Int64, ::Val{:shearxz})::Array{Float64, 3}
    return zeros(Float64, ynum-1, xnum, znum)
end

function grid_array3D(ynum::Int64, xnum::Int64, znum::Int64, ::Val{:shearyz})::Array{Float64, 3}
    return zeros(Float64, ynum, xnum-1, znum)
end

"""
    ScalarArray3DState

Mutable struct representing a 3D scalar array with associated metadata.

# Fields
- `ny::Int64`: Number of nodes in the y-direction
- `nx::Int64`: Number of nodes in the x-direction
- `nz::Int64`: Number of nodes in the z-direction
- `array::Matrix{Float64}`: The 3D array data
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `grid_type::String`: Type of grid ("basic", "pressure", "vx", "vy")
- `description::String`: Description of the array's purpose
- `outform::OutputFormat`: Output formatting specifications
"""
mutable struct ScalarArray3DState <: AbstractEarthBoxArray3D
    ny::Int64
    nx::Int64
    nz::Int64
    array::Array{Float64, 3}
    name::String
    units::String
    grid_type::Symbol
    description::String
    outform::OutputFormat
end

"""
    ScalarArray3DState(number_of_nodes_y, number_of_nodes_x, number_of_nodes_z, name, units, grid_type, description)

Construct a new ScalarArray3DState with the specified parameters.

# Arguments
- `ynum::Int64`: Number of nodes in the y-direction
- `xnum::Int64`: Number of nodes in the x-direction
- `znum::Int64`: Number of nodes in the z-direction
- `name::String`: Name of the array
- `units::String`: Physical units of the array values
- `grid_type::Union{Symbol, String}`: Type of grid (:basic, :pressure, :vx, :vy, :vz, :shearxy, :shearxz, :shearyz)
- `description::String`: Description of the array's purpose

# Returns
- `ScalarArray3DState`: A new instance with initialized array and metadata

# Grid Types
- "basic": Standard grid with dimensions (ny, nx)
- "pressure": Pressure grid with dimensions (ny-1, nx-1, nz-1)
- "vx": Velocity x-component grid with dimensions (ny+1, nx, nz+1)
- "vy": Velocity y-component grid with dimensions (ny, nx+1, nz+1)
- "vz": Velocity z-component grid with dimensions (ny+1, nx+1, nz)
- "shearxy": Shear xy-component grid with dimensions (ny, nx, nz-1)
- "shearxz": Shear xz-component grid with dimensions (ny-1, nx, nz)
- "shearyz": Shear yz-component grid with dimensions (ny, nx-1, nz)
"""
function ScalarArray3DState(
    ynum::Int64,
    xnum::Int64,
    znum::Int64,
    name::String,
    units::String,
    grid_type::Union{Symbol, String},
    description::String
)::ScalarArray3DState
    @assert ynum > 0
    @assert xnum > 0
    @assert znum > 0
    grid_type = Symbol(grid_type)
    array = grid_array3D(ynum, xnum, znum, Val(grid_type))
    outform = OutputFormat(1.0, 0.0, units, name, false, "undefined.dat", name)
    return ScalarArray3DState(
        ynum, xnum, znum, array, name, units, grid_type, description, outform)
end

end # module 