"""
    Functions used to set bc coefficients for prescribed velocity.
"""
module VelocityInflowOutflow

import EarthBox.ModelDataContainer: ModelData
using ...BoundaryVelocity.StickyThickness: get_sticky_thickness
using ...BoundaryVelocity.OutflowGeometry: get_outflow_geometry
using ...BoundaryVelocity.InflowRamp: get_top_and_bottom_of_inflow_ramp, calculate_ramp_factor

"""
    set_inflow_and_outflow_velocity_x_at_left_boundary_ghost_nodes!(
        plate_velocity::Float64,
        inflow_velocity::Float64,
        bc_gridy::Vector{Float64},
        model::ModelData
    )::Nothing

Define coeff. for vx along the left boundary with linear smoothing.

# Arguments
- `sticky_thickness::Float64`: Thickness (m) of the sticky air layer.
- `plate_thickness::Float64`: Thickness (m) of the plate.
- `smoothing_thickness::Float64`: Thickness (m) of the smoothing layer.
- `plate_velocity::Float64`: Velocity (m/s) of the plate.
- `inflow_velocity::Float64`: Velocity (m/s) of the inflow.
- `bc_gridy::Vector{Float64}`: Array of y-coordinates (m) for the staggered grid.

# Updated Array
- `model.bcs.arrays.velocity.bleft.array::Matrix{Float64}`: Array((ynum+1,4))

# Background
The x-component of velocity vx is defined on a staggered grid with ghost nodes.
The following equation is used to define vx at ghost nodes along the left boundary:

vx[i,1] = bleft[i,1] + vx[i,2]*bleft[i,2]

where i is the y-index and the dimensions of vx are (ynum+1, xnum).

Prescribed velocity is implemented by setting vx equal to a user defined value:

vx[i,1] = vx_prescribed

where bleft[i,1] = vx_prescribed and bleft[i,2] = 0. For this case vx_prescribed
is equal to the full_velocity_extension/2.0 until y_linear and then linearly
decreases to zero at the bottom boundary.

# Description of the 2D Staggered Grid Along Left Boundary
```
   j   1            2            3
       vx           vx           vx

i
1 vy     +-----vy-----+-----vy-----+--
         |            |            |
         |            |            |
         vx    p'     vx    p'     vx
         |            |            |
         |            |            |
2 vy     +-----vy-----+-----vy-----+--
         |            |            |
         |            |            |
         vx    p'     vx    p'     vx
         |            |            |
         |            |            |
3 vy     +-----vy-----+-----vy-----+--
         |            |            |
         |            |            |
         vx    p'     vx    p'     vx
         |            |            |
         |            |            |
4 vy     +-----vy-----+-----vy-----+--
         |            |            |
         |            |            |
         vx    p'     vx    p'     vx
         |            |            |
         |            |            |
5 vy     +-----vy-----+-----vy-----+---

       vx           vx           vx
```

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and
continuity equations. Indices i and j apply to basic grid nodes denoted with "+".
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components
of velocity respectively. p' refers to the staggered sub-grid nodes for scaled
pressure with nodes located in the center of basic grid cells. Note that
staggered sub-grids have different indices than the ones displayed in this figure
for the basic grid.
"""
function set_inflow_and_outflow_velocity_x_at_left_boundary_ghost_nodes!(
    plate_velocity::Float64,
    inflow_velocity::Float64,
    bc_gridy::Vector{Float64},
    model::ModelData
)::Nothing
    sticky_thickness_left, _ = get_sticky_thickness(model)
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    y_smooth_top_left, y_smooth_bottom_left = get_top_and_bottom_of_inflow_ramp(
        sticky_thickness_left, plate_thickness, smoothing_thickness)
    ynum = model.grids.parameters.geometry.ynum.value
    bleft = model.bcs.arrays.velocity.bleft.array
    for i in 1:ynum+1
        y_node = bc_gridy[i]
        velocity = calc_velocity_inflow_outflow_along_side(
            y_node, sticky_thickness_left,
            plate_thickness, smoothing_thickness,
            plate_velocity, inflow_velocity,
            y_smooth_top_left, y_smooth_bottom_left
        )
        bleft[i,1] = velocity
        bleft[i,2] = 0.0
    end
    return nothing
end

"""
    set_inflow_and_outflow_velocity_x_at_right_boundary_ghost_nodes!(
        plate_velocity::Float64,
        inflow_velocity::Float64,
        bc_gridy::Vector{Float64},
        model::ModelData
    )::Nothing

Define coeff. for vx along the right boundary with linear smoothing.

# Arguments
- `sticky_thickness::Float64`: Thickness (m) of the sticky air layer.
- `plate_thickness::Float64`: Thickness (m) of the plate.
- `smoothing_thickness::Float64`: Thickness (m) of the smoothing layer.
- `plate_velocity::Float64`: Velocity (m/s) of the plate.
- `inflow_velocity::Float64`: Velocity (m/s) of the inflow.
- `bc_gridy::Vector{Float64}`: Array of y-coordinates (m) for the staggered grid.

# Updated Array
- `model.bcs.arrays.velocity.bright.array::Matrix{Float64}`: Array((ynum+1,4))

# Background
The x-component of velocity vx is defined on a staggered grid with ghost nodes.
The following equation is used to define vx at ghost nodes along the right
boundary:

vx[i,1] = bleft[i,1] + vx[i,2]*bleft[i,2]

where i is the y-index and the dimensions of vx are (ynum+1, xnum).

Prescribed velocity is implemented by setting vx equal to a user defined value:

vx[i,1] = vx_prescribed

where bleft[i,1] = vx_prescribed and bleft[i,2] = 0. For this case vx_prescribed
is equal to the full_velocity_extension/2.0 until y_linear and then linearly
decreases to zero at the bottom boundary.

# Description of the 2D Staggered Grid Along Right Boundary
```
                                j  4           5            6

                                    vx           vx           vx
                                i

                                1---+-----vy-----+-----vy-----+     vy
                                    |            |            |
                                    |            |            |
                                    vx    p'     vx    p'     vx
                                    |            |            |
                                    |            |            |
                                2 ---+-----vy-----+-----vy-----+     vy
                                    |            |            |
                                    |            |            |
                                    vx    p'     vx    p'     vx
                                    |            |            |
                                    |            |            |
                                3 ---+-----vy-----+-----vy-----+     vy
                                    |            |            |
                                    |            |            |
                                    vx     p'     vx    p'     vx
                                    |            |            |
                                    |            |            |
                                4 ---+-----vy-----+-----vy-----+     vy
                                    |            |            |
                                    |            |            |
                                    vx    p'      vx    p'     vx
                                    |            |            |
                                    |            |            |
                                5 ---+-----vy-----+-----vy-----+     vy


                                    vx           vx           vx
```

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and
continuity equations. Indices i and j apply to basic grid nodes denoted with "+".
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components
of velocity respectively. p' refers to the staggered sub-grid nodes for scaled
pressure with nodes located in the center of basic grid cells. Note that
staggered sub-grids have different indices than the ones displayed in this figure
for the basic grid.
"""
function set_inflow_and_outflow_velocity_x_at_right_boundary_ghost_nodes!(
    plate_velocity::Float64,
    inflow_velocity::Float64,
    bc_gridy::Vector{Float64},
    model::ModelData
)::Nothing
    _, sticky_thickness_right = get_sticky_thickness(model)
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    y_smooth_top_right, y_smooth_bottom_right = get_top_and_bottom_of_inflow_ramp(
        sticky_thickness_right, plate_thickness, smoothing_thickness)
    ynum = model.grids.parameters.geometry.ynum.value
    bright = model.bcs.arrays.velocity.bright.array
    for i in 1:ynum+1
        y_node = bc_gridy[i]
        velocity = calc_velocity_inflow_outflow_along_side(
            y_node, sticky_thickness_right,
            plate_thickness, smoothing_thickness,
            plate_velocity, inflow_velocity,
            y_smooth_top_right, y_smooth_bottom_right
        )
        bright[i,1] = velocity
        bright[i,2] = 0.0
    end
    return nothing
end

"""
    calc_velocity_inflow_outflow_along_side(
        y_node::Float64,
        sticky_thickness::Float64,
        plate_thickness::Float64,
        smoothing_thickness::Float64,
        plate_velocity::Float64,
        inflow_velocity::Float64,
        y_inflow_smooth_top::Float64,
        y_inflow_smooth_bottom::Float64
    )::Float64

Calculate velocity for inflow/outflow along side boundaries.

# Arguments
- `y_node::Float64`: Y-coordinate of the node.
- `sticky_thickness::Float64`: Thickness of sticky air layer.
- `plate_thickness::Float64`: Thickness of the plate.
- `smoothing_thickness::Float64`: Thickness of smoothing layer.
- `plate_velocity::Float64`: Velocity of the plate.
- `inflow_velocity::Float64`: Velocity of the inflow.
- `y_inflow_smooth_top::Float64`: Top of inflow smoothing region.
- `y_inflow_smooth_bottom::Float64`: Bottom of inflow smoothing region.

# Returns
- `Float64`: Calculated velocity at the given node.
"""
function calc_velocity_inflow_outflow_along_side(
    y_node::Float64,
    sticky_thickness::Float64,
    plate_thickness::Float64,
    smoothing_thickness::Float64,
    plate_velocity::Float64,
    inflow_velocity::Float64,
    y_inflow_smooth_top::Float64,
    y_inflow_smooth_bottom::Float64
)::Float64
    if y_node > y_inflow_smooth_top
        return inflow_velocity
    elseif y_node < y_inflow_smooth_bottom
        return plate_velocity
    else
        ramp_factor = calculate_ramp_factor(
            y_node, y_inflow_smooth_top, y_inflow_smooth_bottom)
        return inflow_velocity * ramp_factor + plate_velocity * (1.0 - ramp_factor)
    end
end

end # module