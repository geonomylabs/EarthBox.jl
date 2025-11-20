""" Functions used to map markers to basic and pressure grid nodes.

                 j  xn                           xn+1

        i  yn    grid(yn,xn)--------------------grid(yn,xn+1)
                    |           ^                  |
                    |           |                  |
                    |          dy                  |
                    |           |                  |
                    |           v                  |
                    |<----dx--->x marker(ym,xm)    |
                    |                              |
                    |                              |
           yn+1  grid(yn+1,xn)-------------------grid(yn+1, xn+1)

    Figure 1: Schematic showing how surrounding grid nodes are related
    to a marker located within a grid cell. The upper-left grid node
    (yn, xn) is calculated for a given marker location (ym, xm).

"""
module MarkerGridMapping

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs

"""  Find the left-hand index of a point on a 1D grid using bisection.

# Arguments
- `x_marker::Float64': X-location (meters) of marker on 1D grid
- `gridx::Vector{Float64}`: 1D grid of coordinates in meters
- `xnum::Int`: Number of nodes in 1D grid

# Returns
- `ix_upr_left::Int`
    - Index of the left node
"""
@inline function get_index_of_left_node_old(
    x_marker::Float64,
    gridx::Vector{Float64},
    xnum::Int
)::Int32
    # Find horizontal index
    xnmin = 1
    xnmax = xnum
    while (xnmax - xnmin) > 1
        # Subtract 0.5 since ceil(0.5) = 1
        ix_upr_left = ceil(Int32, (xnmax + xnmin)/2 - 0.5)
        #ix_upr_left = (xnmax + xnmin) ÷ 2
        if gridx[ix_upr_left] > x_marker
            xnmax = ix_upr_left
        else
            xnmin = ix_upr_left
        end
    end
    ix_upr_left = xnmin
    # Check horizontal index
    #ix_upr_left = max(ix_upr_left, 1)
    if ix_upr_left < 1
        ix_upr_left = 1
    end
    if ix_upr_left > xnum - 1
        ix_upr_left = xnum - 1
    end
    return ix_upr_left
end

"""  Find the left-hand index of a point on a 1D grid using searchsortedlast.

# Arguments
- `x_marker::Float64': X-location (meters) of marker on 1D grid
- `gridx::Vector{Float64}`: 1D grid of coordinates in meters
- `xnum::Int`: Number of nodes in 1D grid

# Returns
- `ix_upr_left::Int`
    - Index of the left node
"""
@inline function get_index_of_left_node(
    x_marker::Float64,
    gridx::Vector{Float64},
    xnum::Int
)::Int32
    i = searchsortedlast(gridx, x_marker)
    return ifelse(i < 1, 1, ifelse(i > xnum-1, xnum-1, i))
end

"""  Calculate upper-left grid node indices for basic grid.

# Arguments
- `model::ModelData`
    - Model data container

# Updated Arrays
## Updates from group `model.markers.arrays.grid_marker_relationship`:
- `marker_yn.array::Vector{Int32}` (marknum): Vertical index of upper-left basic node
- `marker_xn.array::Vector{Int32}` (marknum): Horizontal index of upper-left basic node
"""
function calc_upper_left_basic_grid_cell_indices_for_markers!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array

    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    marknum = model.markers.parameters.distribution.marknum.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    grid_marker_relationship = model.markers.arrays.grid_marker_relationship
    marker_xn = grid_marker_relationship.marker_xn.array
    marker_yn = grid_marker_relationship.marker_yn.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                marker_xn[imarker] = get_index_of_left_node(x_marker, gridx_b, xnum)
                marker_yn[imarker] = get_index_of_left_node(y_marker, gridy_b, ynum)
            end
        end
    end
    return nothing
end

""" Calculate upper-left grid node indices for vy staggered grid.

# Arguments
- `model::ModelData`
    - Model data container

# Updated Arrays
## Updates from group `model.markers.arrays.grid_marker_relationship`:
- `marker_yn_vy.array::Vector{Int32}` (marknum): Vertical index of upper-left vy staggered grid node
- `marker_xn_vy.array::Vector{Int32}` (marknum): Horizontal index of upper-left vy staggered grid node
"""
function calc_upper_left_vy_grid_cell_indices_for_markers!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    vy_gridx = model.grids.arrays.staggered_vy.gridx_vy.array
    # vy staggered grid has the same y-coordinates as the basic grid
    gridy_b = model.grids.arrays.basic.gridy_b.array

    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    marknum = model.markers.parameters.distribution.marknum.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    grid_marker_relationship = model.markers.arrays.grid_marker_relationship
    marker_xn_vy = grid_marker_relationship.marker_xn_vy.array
    marker_yn_vy = grid_marker_relationship.marker_yn_vy.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                marker_xn_vy[imarker] = get_index_of_left_node(x_marker, vy_gridx, xnum+1)
                marker_yn_vy[imarker] = get_index_of_left_node(y_marker, gridy_b, ynum)
            end
        end
    end
    return nothing
end

"""  Calculate upper-left grid node indices and normalized distances for basic grid.

# Arguments
- `model::ModelData`: Model data container

# Updated Arrays
## Updates from group `model.markers.arrays.grid_marker_relationship`:
- `marker_dx.array::Vector{Float64}` (marknum): Normalized x-distance between marker and upper left grid node
- `marker_dy.array::Vector{Float64}` (marknum): Normalized y-distance between marker and upper left grid node
"""
function calc_normalized_distances_from_upper_left_basic_node!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value

    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array

    xstp_b = model.grids.arrays.basic.xstp_b.array
    ystp_b = model.grids.arrays.basic.ystp_b.array

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array

    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                ix_upper_left = marker_xn[imarker]
                iy_upper_left = marker_yn[imarker]
                marker_dx[imarker] = (x_marker - gridx_b[ix_upper_left])/xstp_b[ix_upper_left]
                marker_dy[imarker] = (y_marker - gridy_b[iy_upper_left])/ystp_b[iy_upper_left]
            end
        end
    end
    return nothing
end

"""  Calculate upper-left grid node indices and normalized distances for vy grid.

# Arguments
- `model::ModelData`: Model data container

# Updated Arrays
## Updates from group `model.markers.arrays.grid_marker_relationship`:
- `marker_dx_vy.array::Vector{Float64}` (marknum): Normalized x-distance between marker and upper left grid node
- `marker_dy_vy.array::Vector{Float64}` (marknum): Normalized y-distance between marker and upper left grid node
"""
function calc_normalized_distances_from_upper_left_vy_node!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value

    vy_gridx = model.grids.arrays.staggered_vy.gridx_vy.array
    gridy_b = model.grids.arrays.basic.gridy_b.array

    xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    ystp_b = model.grids.arrays.basic.ystp_b.array

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    marker_xn_vy = model.markers.arrays.grid_marker_relationship.marker_xn_vy.array
    marker_yn_vy = model.markers.arrays.grid_marker_relationship.marker_yn_vy.array

    marker_dx_vy = model.markers.arrays.grid_marker_relationship.marker_dx_vy.array
    marker_dy_vy = model.markers.arrays.grid_marker_relationship.marker_dy_vy.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                ix_upper_left = marker_xn_vy[imarker]
                iy_upper_left = marker_yn_vy[imarker]
                marker_dx_vy[imarker] = (x_marker - vy_gridx[ix_upper_left])/xstp_vy[ix_upper_left]
                marker_dy_vy[imarker] = (y_marker - gridy_b[iy_upper_left])/ystp_b[iy_upper_left]
            end
        end
    end
    return nothing
end

""" Get upper left node indices from pressure grid for marker.

This functions works for both x- and y-indices.

# Arguments
- `nnodes::Int`:
    - Number of grid nodes in basic grid either for x- or y-direction.
- `index_basic_upr_left::Int`:
    - Index of upper left basic grid node either for x- or y-direction.
- `marker_loc::Float64`:
    - Either x or y location of marker in meters.
- `grid_staggered::Vector{Float64}`:
    - 1D staggered grid array for either x-direction (staggered_vy) or
    y-direction (staggered_vx). Staggered grid arrays have some nodes that
    correspond to pressure grid nodes. For figures illustrating staggered
    grids see the docstring in stokes_continuity_solver.build.stokes_build.

# Returns
- `index_pressure_upr_left::Int`:
    - Index of upper-left node of pressure grid cell.

# Background
Consider the two markers located at m1 and m2 in Figure 1. The
upper-left basic grid node index (xn) in the x-direction for both markers
is 2. If the marker is to the left of staggered grid node 3 (i.e. xn + 1)
then the upper-left index of the pressure grid is 1 or xn-1. However,
if the marker is to the right of staggered grid node 3 then the upper
left pressure grid node index is 2.

      1            2            3            4            5
vy    +-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
1     |     2      |     3      |     4      |     5      |     6
      |            |            |            |            |
      |     P      |     P      |     P      |     P      |
      |     1      |     2      |     3      |     4      |
      |            | m1      m2 |            |            |
vy    +-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy

    Figure 1: Schematic showing basic grid nodes (+), staggered vy grid
    nodes (vy) and pressure grid nodes (P). Numbers correspond to 1D
    x-direction indices of the basic, staggered and pressure grids. Note
    that similar relationships also occurs in the vertical y-direction
    between the staggered-vx grid and the pressure grid.
"""
@inline function upr_left_index_pressure(
    nnodes::Int,
    index_basic_upr_left::Int32,
    marker_loc::Float64,
    grid_staggered::Vector{Float64}
)::Int32
    @inbounds begin
        if marker_loc < grid_staggered[index_basic_upr_left + 1]
            index_pressure_upr_left = index_basic_upr_left - 1
        else
            index_pressure_upr_left = index_basic_upr_left
        end
    end
    index_pressure_upr_left = max(index_pressure_upr_left, 1)
    # Note that the last pressure cell in x-direction has upper left node
    # nnodes-2.
    if index_pressure_upr_left > nnodes-2
        index_pressure_upr_left = nnodes-2
    end
    return index_pressure_upr_left
end

"""  Calculate the normalized distance to upper left node for either x- or y-direction.

# Returns
- `dist_upr_left::Float64`: Distance from marker to upper left node in either 
   x- or y-direction.
"""
@inline function upr_left_dist_pressure(
    index_pressure_upr_left::Int32,
    marker_loc::Float64,
    grid_staggered::Vector{Float64},
    xstp_staggered::Vector{Float64}
)::Float64
    @inbounds dist_upr_left = (
        (marker_loc - grid_staggered[index_pressure_upr_left + 1])
        /xstp_staggered[index_pressure_upr_left + 1]
        )
    return dist_upr_left
end

""" Calculate normalized x-distance to upper left node for basic grid.
"""
@inline function upr_left_x_mapping_basic_grid(
    x_index_upr_left_basic::Int32,
    x_marker::Float64,
    gridx_b::Vector{Float64},
    xstp_b::Vector{Float64}
)::Tuple{Int32,Float64}
    @inbounds dx_upr_left = (x_marker - gridx_b[x_index_upr_left_basic])/xstp_b[x_index_upr_left_basic]
    return x_index_upr_left_basic, dx_upr_left
end

""" Calculate normalized y-distance to upper left node for basic grid.
"""
@inline function upr_left_y_mapping_basic_grid(
    y_index_upr_left_basic::Int32,
    y_marker::Float64,
    gridy_b::Vector{Float64},
    ystp_b::Vector{Float64}
)::Tuple{Int32,Float64}
    @inbounds dy_upr_left = (y_marker - gridy_b[y_index_upr_left_basic])/ystp_b[y_index_upr_left_basic]
    return y_index_upr_left_basic, dy_upr_left
end

""" Calculate normalized x-distance to upper left node for vx grid.
"""
@inline function upr_left_x_mapping_vx_grid(
    x_index_basic_upr_left::Int32,
    x_marker::Float64,
    gridx_b::Vector{Float64},
    xstp_b::Vector{Float64}
)::Tuple{Int32,Float64}
    x_index_upr_left = x_index_basic_upr_left
    @inbounds dx_upr_left = (x_marker - gridx_b[x_index_basic_upr_left])/xstp_b[x_index_basic_upr_left]
    return x_index_upr_left, dx_upr_left
end

""" Calculate normalized y-distance to upper left node for vx grid.
"""
@inline function upr_left_y_mapping_vx_grid(
    y_index_basic_upr_left::Int32,
    y_marker::Float64,
    ynum::Int,
    vxgridy::Vector{Float64},
    ystp_vx::Vector{Float64}
)::Tuple{Int32,Float64}
    y_index_upr_left = y_index_basic_upr_left
    @inbounds begin
        if y_marker > vxgridy[y_index_upr_left+1]
            y_index_upr_left += 1
        end
    end
    if y_index_upr_left > ynum
        y_index_upr_left -= 1
    end
    @inbounds dy_upr_left = (y_marker - vxgridy[y_index_upr_left])/ystp_vx[y_index_upr_left]
    return y_index_upr_left, dy_upr_left
end

""" Calculate normalized x-distance to upper left node for vy grid.
"""
@inline function upr_left_x_mapping_vy_grid(
    x_index_upr_left_basic::Int32,
    x_marker::Float64,
    xnum::Int,
    vygridx::Vector{Float64},
    xstp_vy::Vector{Float64}
)::Tuple{Int32,Float64}
    x_index_upr_left = x_index_upr_left_basic
    if x_marker > vygridx[x_index_upr_left+1]
        x_index_upr_left += 1
    end
    x_index_upr_left = min(x_index_upr_left, xnum)
    # Define and check normalized distances from marker to the upper left VY-node
    @inbounds dx_upr_left = (x_marker - vygridx[x_index_upr_left])/xstp_vy[x_index_upr_left]
    return x_index_upr_left, dx_upr_left
end

""" Calculate normalized y-distance to upper left node for vy grid.
"""
@inline function upr_left_y_mapping_vy_grid(
    y_index_upr_left_basic::Int32,
    y_marker::Float64,
    gridy_b::Vector{Float64},
    ystp_b::Vector{Float64}
)::Tuple{Int32,Float64}
    y_index_upr_left = y_index_upr_left_basic
    @inbounds dy_upr_left = (y_marker - gridy_b[y_index_upr_left_basic])/ystp_b[y_index_upr_left_basic]
    return y_index_upr_left, dy_upr_left
end

end # module 