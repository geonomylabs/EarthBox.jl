module VzGridTools

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function calculate_zgrid_coordinates!(grids3d::Grids3d)::Nothing
    gridz_b = grids3d.arrays.basic.gridz_b.array
    gridz_vz = grids3d.arrays.staggered_vz.gridz_vz.array
    gridz_vz .= gridz_b[1:end]
    return nothing
end

function calculate_xgrid_coordinates!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    gridx_b = grids3d.arrays.basic.gridx_b.array
    xstp_b = grids3d.arrays.basic.xstp_b.array
    gridx_vz = grids3d.arrays.staggered_vz.gridx_vz.array
    vzgrid_x_or_y_coordinates_loop!(xnum, gridx_b, xstp_b, gridx_vz)
    return nothing
end

function calculate_ygrid_coordinates!(grids3d::Grids3d)::Nothing
    ynum = grids3d.parameters.geometry.ynum.value
    gridy_b = grids3d.arrays.basic.gridy_b.array
    ystp_b = grids3d.arrays.basic.ystp_b.array
    gridy_vz = grids3d.arrays.staggered_vz.gridy_vz.array
    vzgrid_x_or_y_coordinates_loop!(ynum, gridy_b, ystp_b, gridy_vz)
    return nothing
end

function vzgrid_x_or_y_coordinates_loop!(
    num::Int64,
    grid_b::Vector{Float64},
    stp_b::Vector{Float64},
    grid_vz::Vector{Float64}
)::Nothing
    grid_vz[1] = grid_b[1] - stp_b[1]/2.0
    grid_vz[num+1] = grid_b[num] + stp_b[num-1]/2.0
    for i in 2:num
        grid_vz[i] = (grid_b[i] + grid_b[i-1])/2.0
    end
    return nothing
end

function calculate_zspacing!(grids3d::Grids3d)::Nothing
    zstp_b = grids3d.arrays.basic.zstp_b.array
    zstp_vz = grids3d.arrays.staggered_vz.zstp_vz.array
    zstp_vz .= zstp_b[1:end]
    return nothing
end

function calculate_xspacing!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    gridx_b = grids3d.arrays.basic.gridx_b.array
    xstp_b = grids3d.arrays.basic.xstp_b.array
    xstp_vz = grids3d.arrays.staggered_vz.xstp_vz.array
    vzgrid_x_or_y_spacing_loop!(xnum, gridx_b, xstp_b, xstp_vz)
    return nothing
end

function calculate_yspacing!(grids3d::Grids3d)::Nothing
    ynum = grids3d.parameters.geometry.ynum.value
    gridy_b = grids3d.arrays.basic.gridy_b.array
    ystp_b = grids3d.arrays.basic.ystp_b.array
    ystp_vz = grids3d.arrays.staggered_vz.ystp_vz.array
    vzgrid_x_or_y_spacing_loop!(ynum, gridy_b, ystp_b, ystp_vz)
    return nothing
end

function vzgrid_x_or_y_spacing_loop!(
    num::Int64,
    grid_b::Vector{Float64},
    stp_b::Vector{Float64},
    stp_vz::Vector{Float64}
)::Nothing
    stp_vz[1] = stp_b[1]
    stp_vz[num] = stp_b[num-1]
    for i in 2:(num-1)
        stp_vz[i] = (grid_b[i+1] - grid_b[i-1])/2.0
    end
    return nothing
end

end # module 