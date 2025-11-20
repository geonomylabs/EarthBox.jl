module GridToMarker

import ..MarkerGridMapping: upr_left_index_pressure
import ..MarkerGridMapping: upr_left_dist_pressure
import ..MarkerGridMapping: get_index_of_left_node
import ..MarkerGridMapping: upr_left_x_mapping_basic_grid
import ..MarkerGridMapping: upr_left_y_mapping_basic_grid

""" Interpolate marker value from basic array without pre-computed mapping.
Grid-marker mapping is calculated by this function.

# Returns
- `marker_value::Float64`: Marker value interpolated from basic grid array.
"""
function get_marker_value_from_basic_grid_without_mapping_input(
    gridx_b::Vector{Float64},
    gridy_b::Vector{Float64},
    xstp_b::Vector{Float64},
    ystp_b::Vector{Float64},
    basic_array::Matrix{Float64},
    x_marker::Float64,
    y_marker::Float64
)::Float64
    xnum = length(gridx_b)
    ynum = length(gridy_b)

    ix_upr_left = get_index_of_left_node(x_marker, gridx_b, xnum)
    iy_upr_left = get_index_of_left_node(y_marker, gridy_b, ynum)

    ix_upr_left, dx_upr_left = upr_left_x_mapping_basic_grid(ix_upr_left, x_marker, gridx_b, xstp_b)
    iy_upr_left, dy_upr_left = upr_left_y_mapping_basic_grid(iy_upr_left, y_marker, gridy_b, ystp_b)

    marker_value = get_marker_value_from_basic_grid_array(
        iy_upr_left,
        ix_upr_left,
        dy_upr_left,
        dx_upr_left,
        basic_array
    )
    return marker_value
end

""" Interpolate marker value from pressure array without pre-computed mapping.
Grid-marker mapping is calculated by this function.

Note that staggered vx and vy grids are used as input.

# Returns
- `marker_value::Float64`: Marker value interpolated from pressure grid array.
"""
function get_marker_value_from_pressure_grid_without_mapping_input(
    gridx_b::Vector{Float64},
    gridy_b::Vector{Float64},
    xstp_b::Vector{Float64},
    ystp_b::Vector{Float64},
    gridy_vx::Vector{Float64},
    ystp_vx::Vector{Float64},
    gridx_vy::Vector{Float64},
    xstp_vy::Vector{Float64},
    pr_array::Matrix{Float64},
    x_marker::Float64,
    y_marker::Float64
)::Float64
    xnum = length(gridx_b)
    ynum = length(gridy_b)

    ix_upr_left = get_index_of_left_node(x_marker, gridx_b, xnum)
    iy_upr_left = get_index_of_left_node(y_marker, gridy_b, ynum)

    ix_upr_left, _dx_upr_left = upr_left_x_mapping_basic_grid(ix_upr_left, x_marker, gridx_b, xstp_b)
    iy_upr_left, _dy_upr_left = upr_left_y_mapping_basic_grid(iy_upr_left, y_marker, gridy_b, ystp_b)

    marker_value = get_marker_value_from_pressure_grid_array(
        ynum, xnum, iy_upr_left, ix_upr_left,
        y_marker, x_marker,
        gridy_vx, ystp_vx, gridx_vy, xstp_vy, pr_array
    )
    return marker_value
end

""" Interpolate marker value from basic grid array.

# Returns
- `marker_value::Float64`: Marker value interpolated from basic grid array.
"""
@inline function get_marker_value_from_basic_grid_array(
    iy_upr_left::Int32,
    ix_upr_left::Int32,
    dy_upr_left::Float64,
    dx_upr_left::Float64,
    basic_array::Matrix{Float64}
)::Float64
    marker_value = get_marker_value(
        iy_upr_left, ix_upr_left, dy_upr_left, dx_upr_left, basic_array
    )
    return marker_value
end

""" Interpolate marker value from pressure grid array.

# Returns
- `marker_value::Float64`: Marker value interpolated from pressure grid array.
"""
@inline function get_marker_value_from_pressure_grid_array(
    ynum::Int,
    xnum::Int,
    iy_upr_left::Int32,
    ix_upr_left::Int32,
    y_marker::Float64,
    x_marker::Float64,
    gridy_vx::Vector{Float64},
    ystp_vx::Vector{Float64},
    gridx_vy::Vector{Float64},
    xstp_vy::Vector{Float64},
    pr_array::Matrix{Float64}
)::Float64
    yn_ul_pr = upr_left_index_pressure(
        ynum,
        iy_upr_left,
        y_marker,
        gridy_vx
    )

    xn_ul_pr = upr_left_index_pressure(
        xnum,
        ix_upr_left,
        x_marker,
        gridx_vy
    )

    dy_ul_pr = upr_left_dist_pressure(
        yn_ul_pr,
        y_marker,
        gridy_vx,
        ystp_vx
    )

    dx_ul_pr = upr_left_dist_pressure(
        xn_ul_pr,
        x_marker,
        gridx_vy,
        xstp_vy
    )

    marker_value = get_marker_value(
        yn_ul_pr,
        xn_ul_pr,
        dy_ul_pr,
        dx_ul_pr,
        pr_array
    )
    return marker_value
end

""" Interpolate grid array values to marker.

A first order bilinear scheme is used for interpolation where markers
are mapped to grid nodes using the upper left node for the grid cell that
contains the marker.

             yn    T[yn][xn]--------------------T[yn][xn+1]
                      |           ^                  |
                      |           |                  |
                      |          dy                  |
                      |           |                  |
                      |           v                  |
                      |<----dx--->o marker           |
                      |                              |
                      |                              |
             yn+1  T[yn+1][xn]-------------------V[yn+1][xn+1]

            Scheme for interpolating grid information from grid array T to
            marker denoted by o.

# Arguments
- `yn_ul::Int`:
    - Y-direction grid node index of the upper-left node of the grid cell
- `xn_ul::Int`:
    - X-direction grid node index of the upper-left node of the grid cell
- `dy_upper_left::Float64`:
    - Y-direction orthogonal distance from the upper left node to the marker
- `dx_upper_left::Float64`:
    - X-direction orthogonal distance from the upper left node to the marker
- `grid_array::Matrix{Float64}`:
    - Grid array to interpolate from

# Returns
- `marker_value::Float64`: Marker value interpolated from grid array
"""
@inline function get_marker_value(
    yn_ul::Int32,
    xn_ul::Int32,
    dy_upper_left::Float64,
    dx_upper_left::Float64,
    grid_array::Matrix{Float64}
)::Float64
    (
        wt_upr_left_node, wt_lwr_left_node, 
        wt_upr_right_node, wt_lwr_right_node
    ) = calc_marker_weights(
        dx_upper_left,
        dy_upper_left
    )
    
    marker_value = 0.0
    @inbounds begin
        marker_value += wt_upr_left_node * grid_array[yn_ul, xn_ul]
        marker_value += wt_lwr_left_node * grid_array[yn_ul+1, xn_ul]
        marker_value += wt_upr_right_node * grid_array[yn_ul, xn_ul+1]
        marker_value += wt_lwr_right_node * grid_array[yn_ul+1, xn_ul+1]
    end
    return marker_value
end

"""  Calculate marker weights for four grid nodes.
"""
@inline function calc_marker_weights(
    dx_upper_left::Float64,
    dy_upper_left::Float64
)::Tuple{Float64,Float64,Float64,Float64}
    wt_upr_left_node = (1.0 - dx_upper_left) * (1.0 - dy_upper_left)
    wt_lwr_left_node = (1.0 - dx_upper_left) * dy_upper_left
    wt_upr_right_node = dx_upper_left * (1.0 - dy_upper_left)
    wt_lwr_right_node = dx_upper_left * dy_upper_left
    return (wt_upr_left_node, wt_lwr_left_node, wt_upr_right_node, wt_lwr_right_node)
end

end # module 