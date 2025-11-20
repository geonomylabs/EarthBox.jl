module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .BasicGridGroup: BasicGrid
import .StaggeredVyGridGroup: StaggeredVyGrid
import .StaggeredVxGridGroup: StaggeredVxGrid
import .StaggeredVzGridGroup: StaggeredVzGrid
import .StaggeredSxyGridGroup: StaggeredSxyGrid
import .StaggeredSxzGridGroup: StaggeredSxzGrid
import .StaggeredSyzGridGroup: StaggeredSyzGrid
import .PressureGridGroup: PressureGrid

mutable struct Arrays <: AbstractArrayCollection
    basic::BasicGrid
    staggered_vy::StaggeredVyGrid
    staggered_vx::StaggeredVxGrid
    staggered_vz::StaggeredVzGrid
    pressure::PressureGrid
    staggered_sxy::StaggeredSxyGrid
    staggered_sxz::StaggeredSxzGrid
    staggered_syz::StaggeredSyzGrid
end

function Arrays(ynum::Int64, xnum::Int64, znum::Int64)::Arrays
    return Arrays(
        BasicGrid(ynum, xnum, znum),
        StaggeredVyGrid(ynum, xnum, znum),
        StaggeredVxGrid(ynum, xnum, znum),
        StaggeredVzGrid(ynum, xnum, znum),
        PressureGrid(ynum, xnum, znum),
        StaggeredSxyGrid(ynum, xnum, znum),
        StaggeredSxzGrid(ynum, xnum, znum),
        StaggeredSyzGrid(ynum, xnum, znum)
    )
end

function get_grid_info(
    data::Arrays,
    grid_type::String
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, String}
    if grid_type == "basic"
        gridx = data.basic.gridx_b.array
        gridy = data.basic.gridy_b.array
        gridz = data.basic.gridz_b.array
        length_units = data.basic.gridx_b.units
    elseif grid_type == "pressure"
        gridx = data.pressure.gridx_pr.array
        gridy = data.pressure.gridy_pr.array
        gridz = data.pressure.gridz_pr.array
        length_units = data.pressure.gridx_pr.units
    elseif grid_type == "staggered_vy"
        gridx = data.staggered_vy.gridx_vy.array
        gridy = data.staggered_vy.gridy_vy.array
        gridz = data.staggered_vy.gridz_vy.array
        length_units = data.staggered_vy.gridx_vy.units
    elseif grid_type == "staggered_vx"
        gridx = data.staggered_vx.gridx_vx.array
        gridy = data.staggered_vx.gridy_vx.array
        gridz = data.staggered_vx.gridz_vx.array
        length_units = data.staggered_vx.gridx_vx.units
    elseif grid_type == "staggered_vz"
        gridx = data.staggered_vz.gridx_vz.array
        gridy = data.staggered_vz.gridy_vz.array
        gridz = data.staggered_vz.gridz_vz.array
        length_units = data.staggered_vz.gridx_vz.units
    elseif grid_type == "staggered_sxy"
        gridx = data.staggered_sxy.gridx_sxy.array
        gridy = data.staggered_sxy.gridy_sxy.array
        gridz = data.staggered_sxy.gridz_sxy.array
        length_units = data.staggered_sxy.gridx_sxy.units
    elseif grid_type == "staggered_sxz"
        gridx = data.staggered_sxz.gridx_sxz.array
        gridy = data.staggered_sxz.gridy_sxz.array
        gridz = data.staggered_sxz.gridz_sxz.array
        length_units = data.staggered_sxz.gridx_sxz.units
    elseif grid_type == "staggered_syz"
        gridx = data.staggered_syz.gridx_syz.array
        gridy = data.staggered_syz.gridy_syz.array
        gridz = data.staggered_syz.gridz_syz.array
        length_units = data.staggered_syz.gridx_syz.units
    else
        throw(ArgumentError("grid_type $grid_type is not valid."))
    end
    return (gridx, gridy, gridz, length_units)
end

end # module 