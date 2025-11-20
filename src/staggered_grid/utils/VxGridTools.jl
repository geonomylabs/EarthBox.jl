module VxGridTools

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function calculate_xgrid_coordinates!(grids3d::Grids3d)::Nothing
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridx_vx = grids3d.arrays.staggered_vx.gridx_vx.array
    gridx_vx .= gridx_b[1:end]
    return nothing
end

function calculate_ygrid_coordinates!(grids::Grids)::Nothing
    ynum = grids.parameters.geometry.ynum.value
    gridy_b = grids.arrays.basic.gridy_b.array
    ystp_b = grids.arrays.basic.ystp_b.array
    gridy_vx = grids.arrays.staggered_vx.gridy_vx.array
    vxgrid_y_or_z_coordinates_loop!(ynum, gridy_b, ystp_b, gridy_vx)
    return nothing
end

function calculate_ygrid_coordinates!(grids3d::Grids3d)::Nothing
    ynum = grids3d.parameters.geometry.ynum.value
    gridy_b = grids3d.arrays.basic.gridy_b.array
    ystp_b = grids3d.arrays.basic.ystp_b.array
    gridy_vx = grids3d.arrays.staggered_vx.gridy_vx.array
    vxgrid_y_or_z_coordinates_loop!(ynum, gridy_b, ystp_b, gridy_vx)
    return nothing
end

function calculate_zgrid_coordinates!(grids3d::Grids3d)::Nothing
    znum = grids3d.parameters.geometry.znum.value
    gridz_b = grids3d.arrays.basic.gridz_b.array
    zstp_b = grids3d.arrays.basic.zstp_b.array
    gridz_vx = grids3d.arrays.staggered_vx.gridz_vx.array
    vxgrid_y_or_z_coordinates_loop!(znum, gridz_b, zstp_b, gridz_vx)
    return nothing
end

function vxgrid_y_or_z_coordinates_loop!(
    num::Int64,
    grid_b::Vector{Float64},
    stp_b::Vector{Float64},
    grid_vx::Vector{Float64}
)::Nothing
    grid_vx[1] = grid_b[1] - stp_b[1]/2.0
    grid_vx[num+1] = grid_b[num] + stp_b[num-1]/2.0
    for i in 2:num
        grid_vx[i] = (grid_b[i] + grid_b[i-1])/2.0
    end
    return nothing
end

function calculate_xspacing!(grids3d::Grids3d)::Nothing
    xstp_b = grids3d.arrays.basic.xstp_b.array
    xstp_vx = grids3d.arrays.staggered_vx.xstp_vx.array
    xstp_vx .= xstp_b[1:end]
    return nothing
end

function calculate_yspacing!(grids::Grids)::Nothing
    ynum = grids.parameters.geometry.ynum.value
    gridy_b = grids.arrays.basic.gridy_b.array
    ystp_b = grids.arrays.basic.ystp_b.array
    ystp_vx = grids.arrays.staggered_vx.ystp_vx.array
    vxgrid_y_or_z_spacing_loop!(ynum, gridy_b, ystp_b, ystp_vx)
    return nothing
end

function calculate_yspacing!(grids3d::Grids3d)::Nothing
    ynum = grids3d.parameters.geometry.ynum.value
    gridy_b = grids3d.arrays.basic.gridy_b.array
    ystp_b = grids3d.arrays.basic.ystp_b.array
    ystp_vx = grids3d.arrays.staggered_vx.ystp_vx.array
    vxgrid_y_or_z_spacing_loop!(ynum, gridy_b, ystp_b, ystp_vx)
    return nothing
end

function calculate_zspacing!(grids3d::Grids3d)::Nothing
    znum = grids3d.parameters.geometry.znum.value
    gridz_b = grids3d.arrays.basic.gridz_b.array
    zstp_b = grids3d.arrays.basic.zstp_b.array
    zstp_vx = grids3d.arrays.staggered_vx.zstp_vx.array
    vxgrid_y_or_z_spacing_loop!(znum, gridz_b, zstp_b, zstp_vx)
    return nothing
end

function vxgrid_y_or_z_spacing_loop!(
    num::Int64,
    grid_b::Vector{Float64},
    stp_b::Vector{Float64},
    stp_vx::Vector{Float64}
)::Nothing
    stp_vx[1] = stp_b[1]
    stp_vx[num] = stp_b[num-1]
    for i in 2:(num-1)
        stp_vx[i] = (grid_b[i+1] - grid_b[i-1])/2.0
    end
    return nothing
end

end # module 