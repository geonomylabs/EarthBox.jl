module NormalizeCoordinates

import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import EarthBox.ModelDataContainer.Grids2dContainer: Grids

function normalize_coordinates_for_basic_node(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d
)::Tuple{Float64, Float64, Float64}
    y = grid.arrays.basic.gridy_b.array[i]
    x = grid.arrays.basic.gridx_b.array[j]
    z = grid.arrays.basic.gridz_b.array[k]
    yy = y / grid.parameters.geometry.ysize.value
    xx = x / grid.parameters.geometry.xsize.value
    zz = z / grid.parameters.geometry.zsize.value
    return yy, xx, zz
end

function normalize_coordinates_for_basic_node(
    i::Int64,
    j::Int64,
    grid::Grids
)::Tuple{Float64, Float64}
    y = grid.arrays.basic.gridy_b.array[i]
    x = grid.arrays.basic.gridx_b.array[j]
    yy = y / grid.parameters.geometry.ysize.value
    xx = x / grid.parameters.geometry.xsize.value
    return yy, xx
end

function normalize_coordinates_for_pressure_node(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d
)::Tuple{Float64, Float64, Float64}
    y = grid.arrays.pressure.gridy_pr.array[i]
    x = grid.arrays.pressure.gridx_pr.array[j]
    z = grid.arrays.pressure.gridz_pr.array[k]
    yy = y / grid.parameters.geometry.ysize.value
    xx = x / grid.parameters.geometry.xsize.value
    zz = z / grid.parameters.geometry.zsize.value
    return yy, xx, zz
end

function normalize_coordinates_for_pressure_node(
    i::Int64,
    j::Int64,
    grid::Grids
)::Tuple{Float64, Float64}
    y = grid.arrays.pressure.gridy_pr.array[i]
    x = grid.arrays.pressure.gridx_pr.array[j]
    yy = y / grid.parameters.geometry.ysize.value
    xx = x / grid.parameters.geometry.xsize.value
    return yy, xx
end

function normalize_coordinates_for_shearxy_node(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d
)::Tuple{Float64, Float64, Float64}
    y = grid.arrays.staggered_sxy.gridy_sxy.array[i]
    x = grid.arrays.staggered_sxy.gridx_sxy.array[j]
    z = grid.arrays.staggered_sxy.gridz_sxy.array[k]
    yy = y / grid.parameters.geometry.ysize.value
    xx = x / grid.parameters.geometry.xsize.value
    zz = z / grid.parameters.geometry.zsize.value
    return yy, xx, zz
end

function normalize_coordinates_for_shearxz_node(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d
)::Tuple{Float64, Float64, Float64}
    y = grid.arrays.staggered_sxz.gridy_sxz.array[i]
    x = grid.arrays.staggered_sxz.gridx_sxz.array[j]
    z = grid.arrays.staggered_sxz.gridz_sxz.array[k]
    yy = y / grid.parameters.geometry.ysize.value
    xx = x / grid.parameters.geometry.xsize.value
    zz = z / grid.parameters.geometry.zsize.value
    return yy, xx, zz
end

function normalize_coordinates_for_shearyz_node(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d
)::Tuple{Float64, Float64, Float64}
    y = grid.arrays.staggered_syz.gridy_syz.array[i]
    x = grid.arrays.staggered_syz.gridx_syz.array[j]
    z = grid.arrays.staggered_syz.gridz_syz.array[k]
    yy = y / grid.parameters.geometry.ysize.value
    xx = x / grid.parameters.geometry.xsize.value
    zz = z / grid.parameters.geometry.zsize.value
    return yy, xx, zz
end

end # module