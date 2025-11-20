module SxyGridTools

import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function calculate_xgrid_coordinates!(grids3d::Grids3d)::Nothing
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridx_sxy = grids3d.arrays.staggered_sxy.gridx_sxy.array
    gridx_sxy .= gridx_b[1:end]
    return nothing
end

function calculate_ygrid_coordinates!(grids3d::Grids3d)::Nothing
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridy_sxy = grids3d.arrays.staggered_sxy.gridy_sxy.array
    gridy_sxy .= gridy_b[1:end]
    return nothing
end

function calculate_zgrid_coordinates!(grids3d::Grids3d)::Nothing
    znum = grids3d.parameters.geometry.znum.value
    gridz_b = grids3d.arrays.basic.gridz_b.array
    zstp_b = grids3d.arrays.basic.zstp_b.array
    gridz_sxy = grids3d.arrays.staggered_sxy.gridz_sxy.array
    sxy_gridz_coordinates_loop!(znum, gridz_b, zstp_b, gridz_sxy)
    return nothing
end

function sxy_gridz_coordinates_loop!(
    znum::Int64,
    gridz_b::Vector{Float64},
    zstp_b::Vector{Float64},
    gridz_sxy::Vector{Float64}
)::Nothing
    gridz_sxy[1] = gridz_b[1] + zstp_b[1]/2.0
    gridz_sxy[end] = gridz_b[end] - zstp_b[end]/2.0
    for i in 2:znum-2
        gridz_sxy[i] = gridz_b[i] + zstp_b[i]/2.0
    end
    return nothing
end

function calculate_xspacing!(grids3d::Grids3d)::Nothing
    xstp_b = grids3d.arrays.basic.xstp_b.array
    xstp_sxy = grids3d.arrays.staggered_sxy.xstp_sxy.array
    xstp_sxy .= xstp_b[1:end]
    return nothing
end

function calculate_yspacing!(grids3d::Grids3d)::Nothing
    ystp_b = grids3d.arrays.basic.ystp_b.array
    ystp_sxy = grids3d.arrays.staggered_sxy.ystp_sxy.array
    ystp_sxy .= ystp_b[1:end]
    return nothing
end

function calculate_zspacing!(grids3d::Grids3d)::Nothing
    znum = grids3d.parameters.geometry.znum.value
    gridz_sxy = grids3d.arrays.staggered_sxy.gridz_sxy.array
    zstp_sxy = grids3d.arrays.staggered_sxy.zstp_sxy.array
    sxy_gridz_spacing_loop!(znum, gridz_sxy, zstp_sxy)
    return nothing
end

function sxy_gridz_spacing_loop!(
    znum::Int64,
    gridz_sxy::Vector{Float64},
    zstp_sxy::Vector{Float64}
)::Nothing
    for i in 1:znum-2
        zstp_sxy[i] = gridz_sxy[i+1] - gridz_sxy[i]
    end
    return nothing
end

end # module 