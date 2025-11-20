module SxzGridTools

import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function calculate_xgrid_coordinates!(grids3d::Grids3d)::Nothing
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridx_sxz = grids3d.arrays.staggered_sxz.gridx_sxz.array
    gridx_sxz .= gridx_b[1:end]
    return nothing
end

function calculate_zgrid_coordinates!(grids3d::Grids3d)::Nothing
    gridz_b = grids3d.arrays.basic.gridz_b.array
    gridz_sxz = grids3d.arrays.staggered_sxz.gridz_sxz.array
    gridz_sxz .= gridz_b[1:end]
    return nothing
end

function calculate_ygrid_coordinates!(grids3d::Grids3d)::Nothing
    ynum = grids3d.parameters.geometry.ynum.value
    gridy_b = grids3d.arrays.basic.gridy_b.array
    ystp_b = grids3d.arrays.basic.ystp_b.array
    gridy_sxz = grids3d.arrays.staggered_sxz.gridy_sxz.array
    sxz_gridy_coordinates_loop!(ynum, gridy_b, ystp_b, gridy_sxz)
    return nothing
end

function sxz_gridy_coordinates_loop!(
    ynum::Int64,
    gridy_b::Vector{Float64},
    ystp_b::Vector{Float64},
    gridy_sxz::Vector{Float64}
)::Nothing
    gridy_sxz[1] = gridy_b[1] + ystp_b[1]/2.0
    gridy_sxz[end] = gridy_b[end] - ystp_b[end]/2.0
    for i in 2:ynum-2
        gridy_sxz[i] = gridy_b[i] + ystp_b[i]/2.0
    end
    return nothing
end

function calculate_xspacing!(grids3d::Grids3d)::Nothing
    xstp_b = grids3d.arrays.basic.xstp_b.array
    xstp_sxz = grids3d.arrays.staggered_sxz.xstp_sxz.array
    xstp_sxz .= xstp_b[1:end]
    return nothing
end

function calculate_zspacing!(grids3d::Grids3d)::Nothing
    zstp_b = grids3d.arrays.basic.zstp_b.array
    zstp_sxz = grids3d.arrays.staggered_sxz.zstp_sxz.array
    zstp_sxz .= zstp_b[1:end]
    return nothing
end

function calculate_yspacing!(grids3d::Grids3d)::Nothing
    ynum = grids3d.parameters.geometry.ynum.value
    gridy_sxz = grids3d.arrays.staggered_sxz.gridy_sxz.array
    ystp_sxz = grids3d.arrays.staggered_sxz.ystp_sxz.array
    sxz_gridy_spacing_loop!(ynum, gridy_sxz, ystp_sxz)
    return nothing
end

function sxz_gridy_spacing_loop!(
    ynum::Int64,
    gridy_sxz::Vector{Float64},
    ystp_sxz::Vector{Float64}
)::Nothing
    for i in 1:ynum-2
        ystp_sxz[i] = gridy_sxz[i+1] - gridy_sxz[i]
    end
    return nothing
end

end # module 