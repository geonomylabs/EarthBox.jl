module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .BasicGridGroup: BasicGrid
import .StaggeredVyGridGroup: StaggeredVyGrid
import .StaggeredVxGridGroup: StaggeredVxGrid
import .PressureGridGroup: PressureGrid

"""
    Arrays <: AbstractArrayCollection

Data structure containing array groups for basic and staggered grids.

# Fields
- `basic::`[`BasicGrid`](@ref): Coordinate and nodal spacing arrays for basic grid
- `staggered_vy::`[`StaggeredVyGrid`](@ref): Coordinate and nodal spacing arrays for staggered Vy grid
- `staggered_vx::`[`StaggeredVxGrid`](@ref): Coordinate and nodal spacing arrays for staggered Vx grid  
- `pressure::`[`PressureGrid`](@ref): Coordinate and nodal spacing arrays for pressure nodes

# Constructor
    Arrays(ynum::Int, xnum::Int)::Arrays

Create a new Arrays collection with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction
"""
mutable struct Arrays <: AbstractArrayCollection
    basic::BasicGrid
    staggered_vy::StaggeredVyGrid
    staggered_vx::StaggeredVxGrid
    pressure::PressureGrid
end

function Arrays(ynum::Int, xnum::Int)::Arrays
    return Arrays(
        BasicGrid(ynum, xnum),
        StaggeredVyGrid(xnum),
        StaggeredVxGrid(ynum),
        PressureGrid(ynum, xnum)
    )
end

function get_grid_info(
    data::Arrays,
    grid_type::String
)::Tuple{Vector{Float64}, Vector{Float64}, String}
    if grid_type == "basic"
        gridx = data.basic.gridx_b.array
        gridy = data.basic.gridy_b.array
        length_units = data.basic.gridx_b.units
    elseif grid_type == "pressure"
        gridx = data.pressure.gridx_pr.array
        gridy = data.pressure.gridy_pr.array
        length_units = data.pressure.gridx_pr.units
    else
        throw(ArgumentError("grid_type $grid_type is not valid."))
    end
    return (gridx, gridy, length_units)
end

end # module 