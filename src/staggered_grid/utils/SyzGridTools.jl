module SyzGridTools

import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function calculate_ygrid_coordinates!(grids3d::Grids3d)::Nothing
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridy_syz = grids3d.arrays.staggered_syz.gridy_syz.array
    gridy_syz .= gridy_b[1:end]
    return nothing
end

function calculate_zgrid_coordinates!(grids3d::Grids3d)::Nothing
    gridz_b = grids3d.arrays.basic.gridz_b.array
    gridz_syz = grids3d.arrays.staggered_syz.gridz_syz.array
    gridz_syz .= gridz_b[1:end]
    return nothing
end

function calculate_xgrid_coordinates!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    gridx_b = grids3d.arrays.basic.gridx_b.array
    xstp_b = grids3d.arrays.basic.xstp_b.array
    gridx_syz = grids3d.arrays.staggered_syz.gridx_syz.array
    syz_gridx_coordinates_loop!(xnum, gridx_b, xstp_b, gridx_syz)
    return nothing
end

function syz_gridx_coordinates_loop!(
    xnum::Int64,
    gridx_b::Vector{Float64},
    xstp_b::Vector{Float64},
    gridx_syz::Vector{Float64}
)::Nothing
    gridx_syz[1] = gridx_b[1] + xstp_b[1]/2.0
    gridx_syz[end] = gridx_b[end] - xstp_b[end]/2.0
    for i in 2:xnum-2
        gridx_syz[i] = gridx_b[i] + xstp_b[i]/2.0
    end
    return nothing
end

function calculate_yspacing!(grids3d::Grids3d)::Nothing
    ystp_b = grids3d.arrays.basic.ystp_b.array
    ystp_syz = grids3d.arrays.staggered_syz.ystp_syz.array
    ystp_syz .= ystp_b[1:end]
    return nothing
end

function calculate_zspacing!(grids3d::Grids3d)::Nothing
    zstp_b = grids3d.arrays.basic.zstp_b.array
    zstp_syz = grids3d.arrays.staggered_syz.zstp_syz.array
    zstp_syz .= zstp_b[1:end]
    return nothing
end

function calculate_xspacing!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    gridx_syz = grids3d.arrays.staggered_syz.gridx_syz.array
    xstp_syz = grids3d.arrays.staggered_syz.xstp_syz.array
    syz_gridx_spacing_loop!(xnum, gridx_syz, xstp_syz)
    return nothing
end

function syz_gridx_spacing_loop!(
    xnum::Int64,
    gridx_syz::Vector{Float64},
    xstp_syz::Vector{Float64}
)::Nothing
    for i in 1:xnum-2
        xstp_syz[i] = gridx_syz[i+1] - gridx_syz[i]
    end
    return nothing
end

end # module 