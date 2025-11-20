module Spacing

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function update_average_spacing_of_basic_grid!(grids::Grids)::Nothing
    xnum = grids.parameters.geometry.xnum.value
    ynum = grids.parameters.geometry.ynum.value
    xsize = grids.parameters.geometry.xsize.value
    ysize = grids.parameters.geometry.ysize.value

    xstpavg = xsize / (xnum - 1)
    ystpavg = ysize / (ynum - 1)

    grids.parameters.geometry.xstpavg.value = xstpavg
    grids.parameters.geometry.ystpavg.value = ystpavg
    return nothing
end

function update_average_spacing_of_basic_grid!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    ynum = grids3d.parameters.geometry.ynum.value
    znum = grids3d.parameters.geometry.znum.value
    xsize = grids3d.parameters.geometry.xsize.value
    ysize = grids3d.parameters.geometry.ysize.value
    zsize = grids3d.parameters.geometry.zsize.value

    xstpavg = xsize / (xnum - 1)
    ystpavg = ysize / (ynum - 1)
    zstpavg = zsize / (znum - 1)

    grids3d.parameters.geometry.xstpavg.value = xstpavg
    grids3d.parameters.geometry.ystpavg.value = ystpavg
    grids3d.parameters.geometry.zstpavg.value = zstpavg
    return nothing
end

function calculate_average_grid_spacing_using_grid_arrays(
    grids::Grids
)::Tuple{Float64, Float64}
    xnum = grids.parameters.geometry.xnum.value
    ynum = grids.parameters.geometry.ynum.value
    gridx_b = grids.arrays.basic.gridx_b.array
    gridy_b = grids.arrays.basic.gridy_b.array

    xstpavg = (gridx_b[xnum] - gridx_b[1]) / (xnum - 1)
    ystpavg = (gridy_b[ynum] - gridy_b[1]) / (ynum - 1)
    return xstpavg, ystpavg
end

function calculate_average_grid_spacing_using_grid_arrays(
    grids3d::Grids3d
)::Tuple{Float64, Float64, Float64}
    xnum = grids3d.parameters.geometry.xnum.value
    ynum = grids3d.parameters.geometry.ynum.value
    znum = grids3d.parameters.geometry.znum.value
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridz_b = grids3d.arrays.basic.gridz_b.array

    xstpavg = (gridx_b[xnum] - gridx_b[1]) / (xnum - 1)
    ystpavg = (gridy_b[ynum] - gridy_b[1]) / (ynum - 1)
    zstpavg = (gridz_b[znum] - gridz_b[1]) / (znum - 1)
    return xstpavg, ystpavg, zstpavg
end

function update_spacing_of_basic_grid!(grids::Grids)::Nothing
    xnum = grids.parameters.geometry.xnum.value
    ynum = grids.parameters.geometry.ynum.value
    gridx_b = grids.arrays.basic.gridx_b.array
    gridy_b = grids.arrays.basic.gridy_b.array
    xstp_b = grids.arrays.basic.xstp_b.array
    ystp_b = grids.arrays.basic.ystp_b.array
    
    update_basic_1Dgrid_spacing!(xnum, gridx_b, xstp_b)
    update_basic_1Dgrid_spacing!(ynum, gridy_b, ystp_b)
    return nothing
end

function update_spacing_of_basic_grid!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    ynum = grids3d.parameters.geometry.ynum.value
    znum = grids3d.parameters.geometry.znum.value
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridz_b = grids3d.arrays.basic.gridz_b.array
    xstp_b = grids3d.arrays.basic.xstp_b.array
    ystp_b = grids3d.arrays.basic.ystp_b.array
    zstp_b = grids3d.arrays.basic.zstp_b.array
    
    update_basic_1Dgrid_spacing!(xnum, gridx_b, xstp_b)
    update_basic_1Dgrid_spacing!(ynum, gridy_b, ystp_b)
    update_basic_1Dgrid_spacing!(znum, gridz_b, zstp_b)
    return nothing
end

function get_basic_grid_spacing_arrays!(
    gridx_b::Vector{Float64},
    gridy_b::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}}
    xstp, ystp = initialize_spacing_arrays(gridx_b, gridy_b)
    xnum = length(gridx_b)
    ynum = length(gridy_b)
    update_basic_1Dgrid_spacing!(xnum, gridx_b, xstp)
    update_basic_1Dgrid_spacing!(ynum, gridy_b, ystp)
    return xstp, ystp
end

function get_basic_3dgrid_spacing_arrays!(
    gridx_b::Vector{Float64},
    gridy_b::Vector{Float64},
    gridz_b::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    xstp_b, ystp_b, zstp_b = initialize_spacing_arrays(gridx_b, gridy_b, gridz_b)
    xnum = length(gridx_b)
    ynum = length(gridy_b)
    znum = length(gridz_b)
    update_basic_1Dgrid_spacing!(xnum, gridx_b, xstp_b)
    update_basic_1Dgrid_spacing!(ynum, gridy_b, ystp_b)
    update_basic_1Dgrid_spacing!(znum, gridz_b, zstp_b)
    return xstp_b, ystp_b, zstp_b
end

function update_basic_1Dgrid_spacing!(
    num::Int,
    grid::Vector{Float64},
    stp1::Vector{Float64}
)::Nothing
    for i in 1:(num-1)
        stp1[i] = grid[i+1] - grid[i]
    end
    return nothing
end

function initialize_spacing_arrays(
    gridx_b::Vector{Float64},
    gridy_b::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}}
    xnum = length(gridx_b)
    ynum = length(gridy_b)
    xstp_b = zeros(Float64, xnum - 1)
    ystp_b = zeros(Float64, ynum - 1)
    return xstp_b, ystp_b
end

function initialize_spacing_arrays(
    gridx_b::Vector{Float64},
    gridy_b::Vector{Float64},
    gridz_b::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    xnum = length(gridx_b)
    ynum = length(gridy_b)
    znum = length(gridz_b)
    xstp_b = zeros(Float64, xnum - 1)
    ystp_b = zeros(Float64, ynum - 1)
    zstp_b = zeros(Float64, znum - 1)
    return xstp_b, ystp_b, zstp_b
end

end # module Spacing

