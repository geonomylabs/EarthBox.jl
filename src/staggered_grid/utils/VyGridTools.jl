module VyGridTools

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function calculate_ygrid_coordinates!(grids3d::Grids3d)::Nothing
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridy_vy = grids3d.arrays.staggered_vy.gridy_vy.array
    gridy_vy .= gridy_b[1:end]
    return nothing
end

function calculate_xgrid_coordinates!(grids::Grids)::Nothing
    xnum = grids.parameters.geometry.xnum.value
    gridx_b = grids.arrays.basic.gridx_b.array
    xstp_b = grids.arrays.basic.xstp_b.array
    gridx_vy = grids.arrays.staggered_vy.gridx_vy.array
    vygrid_x_or_z_coordinates_loop!(xnum, gridx_b, xstp_b, gridx_vy)
    return nothing
end

function calculate_xgrid_coordinates!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    gridx_b = grids3d.arrays.basic.gridx_b.array
    xstp_b = grids3d.arrays.basic.xstp_b.array
    gridx_vy = grids3d.arrays.staggered_vy.gridx_vy.array
    vygrid_x_or_z_coordinates_loop!(xnum, gridx_b, xstp_b, gridx_vy)
    return nothing
end

function calculate_zgrid_coordinates!(grids3d::Grids3d)::Nothing
    znum = grids3d.parameters.geometry.znum.value
    gridz_b = grids3d.arrays.basic.gridz_b.array
    zstp_b = grids3d.arrays.basic.zstp_b.array
    gridz_vy = grids3d.arrays.staggered_vy.gridz_vy.array
    vygrid_x_or_z_coordinates_loop!(znum, gridz_b, zstp_b, gridz_vy)
    return nothing
end

function vygrid_x_or_z_coordinates_loop!(
    num::Int64,
    grid_b::Vector{Float64},
    stp_b::Vector{Float64},
    grid_vy::Vector{Float64}
)::Nothing
    grid_vy[1] = grid_b[1] - stp_b[1]/2.0
    grid_vy[num+1] = grid_b[num] + stp_b[num-1]/2.0
    for i in 2:num
        grid_vy[i] = (grid_b[i] + grid_b[i-1])/2.0
    end
    return nothing
end

function calculate_yspacing!(grids3d::Grids3d)::Nothing
    ystp_b = grids3d.arrays.basic.ystp_b.array
    ystp_vy = grids3d.arrays.staggered_vy.ystp_vy.array
    ystp_vy .= ystp_b[1:end]
    return nothing
end

function calculate_xspacing!(grids::Grids)::Nothing
    xnum = grids.parameters.geometry.xnum.value
    gridx_b = grids.arrays.basic.gridx_b.array
    xstp_b = grids.arrays.basic.xstp_b.array
    xstp_vy = grids.arrays.staggered_vy.xstp_vy.array
    vygrid_x_or_z_spacing_loop!(xnum, gridx_b, xstp_b, xstp_vy)
    return nothing
end

function calculate_xspacing!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    gridx_b = grids3d.arrays.basic.gridx_b.array
    xstp_b = grids3d.arrays.basic.xstp_b.array
    xstp_vy = grids3d.arrays.staggered_vy.xstp_vy.array
    vygrid_x_or_z_spacing_loop!(xnum, gridx_b, xstp_b, xstp_vy)
    return nothing
end

function calculate_zspacing!(grids3d::Grids3d)::Nothing
    znum = grids3d.parameters.geometry.znum.value
    gridz_b = grids3d.arrays.basic.gridz_b.array
    zstp_b = grids3d.arrays.basic.zstp_b.array
    zstp_vy = grids3d.arrays.staggered_vy.zstp_vy.array
    vygrid_x_or_z_spacing_loop!(znum, gridz_b, zstp_b, zstp_vy)
    return nothing
end

function vygrid_x_or_z_spacing_loop!(
    num::Int64,
    grid_b::Vector{Float64},
    stp_b::Vector{Float64},
    stp_vy::Vector{Float64}
)::Nothing
    stp_vy[1] = stp_b[1]
    stp_vy[num] = stp_b[num-1]
    for i in 2:(num-1)
        stp_vy[i] = (grid_b[i+1] - grid_b[i-1])/2.0
    end
    return nothing
end

end # module 