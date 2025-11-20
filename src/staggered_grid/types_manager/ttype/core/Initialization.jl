module Initialization

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import ....Dimensions: update_grid_dimensions!
import ....Spacing: update_average_spacing_of_basic_grid!

function initialize!(grids::Grids)
    update_average_spacing_of_basic_grid!(grids)
    calculate_gridx_coordinates!(grids)
    calculate_gridy_coordinates!(grids)
    update_grid_dimensions!(grids)
end

function calculate_gridx_coordinates!(grids::Grids)
    # Initial x-grid location of high-resolution area
    xo_highres = grids.parameters.refinement.xo_highres.value
    # Final x-grid location of high-resolution area
    xf_highres = grids.parameters.refinement.xf_highres.value
    # X-direction spacing in high-resolution area
    dx_highres = grids.parameters.refinement.dx_highres.value
    # Number of nodes of basic grid in x-direction
    xnum = grids.parameters.geometry.xnum.value
    # Size of model in x-direction
    xsize = grids.parameters.geometry.xsize.value
    # X-coordinates of basic grid in x-direction
    gridx_b = grids.arrays.basic.gridx_b.array

    ixo_highres = get_left_edge_x_index_of_highres_zone(
        xnum, xo_highres, xf_highres, dx_highres)
    ixf_highres = calculate_final_xgrid_index_of_highres(
        ixo_highres, xo_highres, xf_highres, dx_highres)

    # Define position of first and last nodal points
    gridx_b[1] = 0.0
    gridx_b[xnum] = xsize

    internal_high_res_xgrid!(
        ixo_highres, ixf_highres, xo_highres, dx_highres, gridx_b)

    # Distance to be covered by right-side non-uniform grid
    dist = xsize - gridx_b[ixf_highres]
    define_xgrid_to_the_right_of_high_res_region!(
        xnum, ixf_highres, dist, dx_highres, gridx_b)

    # Distance to be covered by left-side non-uniform grid
    dist = gridx_b[ixo_highres] - gridx_b[1]
    define_xgrid_to_the_left_of_high_res_region!(
        ixo_highres, dist, dx_highres, gridx_b)
end

function get_left_edge_x_index_of_highres_zone(
    xnum::Int64,
    xo_highres::Float64,
    xf_highres::Float64,
    dx_highres::Float64
)::Int64
    nnodes_hr = (xf_highres - xo_highres)/dx_highres + 1
    ixo_highres = Int64(floor((xnum - nnodes_hr)/2.0)) + 1 # +1 to account for 1-based indexing
    return ixo_highres
end

function calculate_final_xgrid_index_of_highres(
    ixo_highres::Int64,
    xo_highres::Float64,
    xf_highres::Float64,
    dx_highres::Float64
)::Int64
    ixf_highres = ixo_highres + Int64(floor((xf_highres - xo_highres)/dx_highres))
    return ixf_highres
end

function define_xgrid_to_the_left_of_high_res_region!(
    ixo_highres::Int64,
    dist::Float64,
    dx_highres::Float64,
    gridx_b::Vector{Float64}
)
    nsteps = ixo_highres - 1
    smoothing_factor = calculate_grid_smoothing_factor(dist, dx_highres, nsteps)
    for i in 2:nsteps
        gridx_b[i] = gridx_b[i-1] + dx_highres * smoothing_factor^(convert(Float64, nsteps + 2 - i))
    end
end

function define_xgrid_to_the_right_of_high_res_region!(
    xnum::Int64,
    ixf_highres::Int64,
    dist::Float64,
    dx_highres::Float64,
    gridx_b::Vector{Float64}
)
    nsteps = xnum - ixf_highres
    smoothing_factor = calculate_grid_smoothing_factor(dist, dx_highres, nsteps)
    for i in (ixf_highres+1):xnum
        gridx_b[i] = gridx_b[i-1] + dx_highres * smoothing_factor^(convert(Float64, i - ixf_highres))
    end
end

function internal_high_res_xgrid!(
    ixo_highres::Int64,
    ixf_highres::Int64,
    xo_highres::Float64,
    dx_highres::Float64,
    gridx_b::Vector{Float64}
)
    gridx_b[ixo_highres] = xo_highres
    for i in (ixo_highres+1):ixf_highres
        gridx_b[i] = gridx_b[i-1] + dx_highres
    end
end

function calculate_gridy_coordinates!(grids::Grids)
    yf_highres = grids.parameters.refinement.yf_highres.value
    dy_highres = grids.parameters.refinement.dy_highres.value
    ynum = grids.parameters.geometry.ynum.value
    ysize = grids.parameters.geometry.ysize.value
    gridy_b = grids.arrays.basic.gridy_b.array

    gridy_b[ynum] = ysize
    iyf_highres = Int64(floor(yf_highres/dy_highres)) + 1

    internal_high_res_ygrid!(iyf_highres, dy_highres, gridy_b)

    # Distance to be covered by non-uniform grid
    dist = ysize - gridy_b[iyf_highres]
    define_ygrid_below_high_res_region!(iyf_highres, ynum, dy_highres, dist, gridy_b)
end

function internal_high_res_ygrid!(
    iyf_highres::Int64,
    dy_highres::Float64,
    gridy_b::Vector{Float64}
)
    for i in 2:iyf_highres
        gridy_b[i] = gridy_b[i-1] + dy_highres
    end
end

function define_ygrid_below_high_res_region!(
    iyf_highres::Int64,
    ynum::Int64,
    dy_highres::Float64,
    dist::Float64,
    gridy_b::Vector{Float64}
)
    nsteps = ynum - iyf_highres
    smoothing_factor = calculate_grid_smoothing_factor(dist, dy_highres, nsteps)

    for i in (iyf_highres+1):ynum
        gridy_b[i] = gridy_b[i-1] + dy_highres * smoothing_factor^(convert(Float64, i - iyf_highres))
    end
end

function calculate_grid_smoothing_factor(
    dist::Float64, dw_highres::Float64,
    nsteps::Int64
)::Float64
    smoothing_factor = 1.1
    for _ in 1:100
        smoothing_factor = (
            1.0 
            + dist/dw_highres*(1.0 - 1.0/smoothing_factor)
            )^(1.0/convert(Float64, nsteps))
    end
    return smoothing_factor
end

end # module