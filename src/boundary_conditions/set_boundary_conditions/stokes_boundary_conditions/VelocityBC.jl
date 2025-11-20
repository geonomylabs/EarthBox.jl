module VelocityBC

import EarthBox.ModelDataContainer: ModelData

"""
    set_velocity_y_at_top_boundary_ghost_nodes!(vy::Float64, model::ModelData)::Nothing

Define coefficients for prescribing vy along the top boundary.

# Updated arrays from group bcs.arrays.velocity
- btop.array::Matrix{Float64}

# Background
The y-component of velocity vy is defined on a staggered grid with ghost nodes. The 
following equation is used to define vy ghost nodes along the top boundary:

    vy[1,j] = btop[j,3] + vy[2,j]*btop[j,4]

where j is the x-index and the dimensions of vy are (ynum, xnum+1).

Prescribed velocity is implemented by setting vy at ghost node equal to a user 
defined value:

    vy[1,j] = vy_prescribed

where btop[j,3] = vy_prescribed and btop[j,4] = 0.

# Description of the 2D Staggered Grid Along Top Boundary
       j   1            2            3             4           5            6
           vx           vx           vx           vx           vx           vx
  i
  1 vy     +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
           |            |            |            |            |            |
           |            |            |            |            |            |
           vx    p'     vx    p'     vx    p'     vx    p'     vx    p'     vx
           |            |            |            |            |            |
           |            |            |            |            |            |
  2 vy     +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
           |            |            |            |            |            |

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and 
continuity equations. Indices i and j apply to basic grid nodes denoted with "+". 
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components of 
velocity respectively. p' refers to the staggered sub-grid nodes for scaled pressure 
with nodes located in the center of basic grid cells. Note that staggered sub-grids 
have different indices than the ones displayed in this figure for the basic grid.
"""
function set_velocity_y_at_top_boundary_ghost_nodes!(vy::Float64, model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    btop = model.bcs.arrays.velocity.btop.array
    for j in 1:xnum+1
        btop[j,3] = vy
        btop[j,4] = 0.0
    end
    return nothing
end

"""
    set_velocity_y_at_bottom_boundary_ghost_nodes!(vy::Float64, model::ModelData)::Nothing

Define coefficients for prescribing vy along the bottom boundary.

# Updated arrays from group bcs.arrays.velocity
- bbottom.array::Matrix{Float64}

# Background
The y-component of velocity vy is defined on a staggered grid with ghost nodes. The 
following equation is used to define vy ghost nodes located along the bottom 
boundary:

    vy[ynum,j] = bbottom[j,3] + vy[ynum-1,j]*bbottom[j,4]

where j is the x-index and the dimensions of vy are (ynum, xnum+1).

Prescribed velocity is implemented by setting vy at ghost node equal to a user 
defined value:

    vy[ynum,j] = vy_prescribed

where bbottom[j,3] = vy_prescribed and bbottom[j,4] = 0.

# Description of the 2D Staggered Grid Along Bottom Boundary
       j   1            2            3             4           5            6
  i
           |            |            |            |            |            |
  3 vy     +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
           |            |            |            |            |            |
           |            |            |            |            |            |
           vx    p'     vx    p'     vx    p'     vx    p'     vx    p'     vx
           |            |            |            |            |            |
           |            |            |            |            |            |
  4 vy     +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
           |            |            |            |            |            |
           |            |            |            |            |            |
           vx    p'     vx    p'     vx    p'     vx    p'     vx    p'     vx
           |            |            |            |            |            |
           |            |            |            |            |            |
  5 vy     +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy

           vx           vx           vx           vx           vx           vx

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and 
continuity equations. Indices i and j apply to basic grid nodes denoted with "+". 
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components of 
velocity respectively. p' refers to the staggered sub-grid nodes for scaled pressure 
with nodes located in the center of basic grid cells. Note that staggered sub-grids 
have different indices than the ones displayed in this figure for the basic grid.
"""
function set_velocity_y_at_bottom_boundary_ghost_nodes!(vy::Float64, model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    bbottom = model.bcs.arrays.velocity.bbottom.array
    for j in 1:xnum+1
        bbottom[j,3] = vy
        bbottom[j,4] = 0.0
    end
    return nothing
end

"""
    set_velocity_x_at_left_boundary_ghost_nodes!(vx::Float64, model::ModelData)::Nothing

Define coefficients for prescribing vx along the left boundary.

# Updated arrays from group bcs.arrays.velocity
- bleft.array::Matrix{Float64}

# Background
The x-component of velocity vx is defined on a staggered grid with ghost nodes. The 
following equation is used to define vx at ghost nodes along the left boundary:

    vx[i,1] = bleft[i,1] + vx[i,2]*bleft[i,2]

where i is the y-index and the dimensions of vx are (ynum+1, xnum).

Prescribed velocity is implemented by setting vx equal to a user defined value:

    vx[i,1] = vx_prescribed

where bleft[i,1] = vx_prescribed and bleft[i,2] = 0.

# Description of the 2D Staggered Grid Along Left Boundary
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
  4 vy     +-----vy-----+-----vy-----+---
           |            |            |
           |            |            |
           vx    p'     vx    p'     vx
           |            |            |
           |            |            |
  5 vy     +-----vy-----+-----vy-----+---

           vx           vx           vx

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and 
continuity equations. Indices i and j apply to basic grid nodes denoted with "+". 
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components of 
velocity respectively. p' refers to the staggered sub-grid nodes for scaled pressure 
with nodes located in the center of basic grid cells. Note that staggered sub-grids 
have different indices than the ones displayed in this figure for the basic grid.
"""
function set_velocity_x_at_left_boundary_ghost_nodes!(
    vx::Float64,
    model::ModelData
)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bleft = model.bcs.arrays.velocity.bleft.array
    for i in 1:ynum+1
        bleft[i,1] = vx
        bleft[i,2] = 0.0
    end
    return nothing
end

"""
    set_linear_velocity_x_at_left_boundary_ghost_nodes!(
        y_linear::Float64,
        velocity_extension::Float64,
        bc_gridy::Vector{Float64},
        model::ModelData
    )::Nothing

Define coefficients for vx along the left boundary with linear smoothing.

# Arguments
- `y_linear`: Depth (m) at which a linear decrease in velocity is applied.
- `bc_gridy`: Array of y-coordinates (m) for the staggered grid.
- `velocity_extension`: Extension of the velocity (m/s).

# Updated arrays from group bcs.arrays.velocity
- bleft.array::Matrix{Float64}

# Background
The x-component of velocity vx is defined on a staggered grid with ghost nodes. The 
following equation is used to define vx at ghost nodes along the left boundary:

    vx[i,1] = bleft[i,1] + vx[i,2]*bleft[i,2]

where i is the y-index and the dimensions of vx are (ynum+1, xnum).

Prescribed velocity is implemented by setting vx equal to a user defined value:

    vx[i,1] = vx_prescribed

where bleft[i,1] = vx_prescribed and bleft[i,2] = 0. For this case vx_prescribed is 
equal to the full_velocity_extension/2.0 until y_linear and then linearly decreases 
to zero at the bottom boundary.
"""
function set_linear_velocity_x_at_left_boundary_ghost_nodes!(
    y_linear::Float64,
    velocity_extension::Float64,
    bc_gridy::Vector{Float64},
    model::ModelData
)::Nothing
    ysize = model.grids.parameters.geometry.ysize.value
    ynum = model.grids.parameters.geometry.ynum.value
    bleft = model.bcs.arrays.velocity.bleft.array
    for i in 1:ynum+1
        yy = bc_gridy[i]
        velbc = calc_velocity_depth_dependent(yy, y_linear, ysize, velocity_extension)
        bleft[i,1] = -velbc
        bleft[i,2] = 0.0
    end
    return nothing
end

"""
    set_linear_velocity_x_at_right_boundary_ghost_nodes!(
        y_linear::Float64,
        velocity_extension::Float64,
        bc_gridy::Vector{Float64},
        model::ModelData
    )::Nothing

Define coefficients for vx along the right boundary with linear smoothing.

# Arguments
- `y_linear`: Depth (m) at which a linear decrease in velocity is applied.
- `bc_gridy`: Array of y-coordinates (m) for the staggered grid.
- `velocity_extension`: Extension velocity (m/s).

# Updated arrays from group bcs.arrays.velocity
- bright.array::Matrix{Float64}

# Background
The x-component of velocity vx is defined on a staggered grid with ghost nodes. The 
following equation is used to define vx at ghost nodes along the right boundary:

    vx[i,1] = bleft[i,1] + vx[i,2]*bleft[i,2]

where i is the y-index and the dimensions of vx are (ynum+1, xnum).

Prescribed velocity is implemented by setting vx equal to a user defined value:

    vx[i,1] = vx_prescribed

where bleft[i,1] = vx_prescribed and bleft[i,2] = 0. For this case vx_prescribed is 
equal to the full_velocity_extension/2.0 until y_linear and then linearly decreases 
to zero at the bottom boundary.
"""
function set_linear_velocity_x_at_right_boundary_ghost_nodes!(
    y_linear::Float64,
    velocity_extension::Float64,
    bc_gridy::Vector{Float64},
    model::ModelData
)::Nothing
    ysize = model.grids.parameters.geometry.ysize.value
    ynum = model.grids.parameters.geometry.ynum.value
    bright = model.bcs.arrays.velocity.bright.array
    for i in 1:ynum+1
        yy = bc_gridy[i]
        velbc = calc_velocity_depth_dependent(yy, y_linear, ysize, velocity_extension)
        bright[i,1] = velbc
        bright[i,2] = 0.0
    end
    return nothing
end

"""
    calc_velocity_depth_dependent(
        yy::Float64,
        y_linear::Float64,
        ysize::Float64,
        velocity_extension::Float64
    )::Float64

Calculate depth-dependent velocity along side boundary.

# Arguments
- `yy`: Depth (m) at which velocity is calculated.
- `y_linear`: Depth (m) at which a linear decrease in velocity is applied.
- `ysize`: Depth (m) of the model.
- `velocity_extension`: Extension velocity (m/s).
"""
function calc_velocity_depth_dependent(
    yy::Float64,
    y_linear::Float64,
    ysize::Float64,
    velocity_extension::Float64
)::Float64
    velocity_extension = abs(velocity_extension)
    if yy <= y_linear
        velbc = velocity_extension
    else
        zz = yy - y_linear
        velbc = velocity_extension - velocity_extension/(ysize - y_linear)*zz
    end
    return velbc
end

""" Define coefficients for prescribing vx along the right boundary.

# Updated Array
model.bcs.arrays.velocity.bright.array: Array((ynum+1,4), dtype=np.float64)


# Background
The x-component of velocity vx is defined on a staggered grid with ghost
nodes. The following equation is used to define vx at ghost nodes along the
right boundary:

vx[i,xnum] = bright[i,0] + vx[i,xnum-1]*bright[i,1]

where i is the y-index and the dimensions of vx are (ynum+1, xnum).

Prescribed velocity is implemented by setting vx equal to a user defined
value:

vx[i,xnum] = vx_prescribed

where i is the y-index and vx_prescribed is equal to bright[i,0].

# Description of the 2D Staggered Grid Along Right Boundary

                            j  3           4            5

                                vx           vx           vx
                            i

                            0---+-----vy-----+-----vy-----+     vy
                                |            |            |
                                |            |            |
                                vx    p'     vx    p'     vx
                                |            |            |
                                |            |            |
                            1 ---+-----vy-----+-----vy-----+     vy
                                |            |            |
                                |            |            |
                                vx    p'     vx    p'     vx
                                |            |            |
                                |            |            |
                            2 ---+-----vy-----+-----vy-----+     vy
                                |            |            |
                                |            |            |
                                vx     p'     vx    p'     vx
                                |            |            |
                                |            |            |
                            3 ---+-----vy-----+-----vy-----+     vy
                                |            |            |
                                |            |            |
                                vx    p'      vx    p'     vx
                                |            |            |
                                |            |            |
                            4 ---+-----vy-----+-----vy-----+     vy


                                vx           vx           vx

    Figure 1: 2D staggered grid and numbering scheme used for discretizing
    Stokes and continuity equations. Indices i and j apply to basic grid
    nodes denoted with "+". Symbols vx and vy refer to staggered sub-grid
    nodes for the x- and y-components of velocity respectively. p' refers
    to the staggered sub-grid nodes for scaled pressure with nodes located
    in the center of basic grid cells. Note that staggered sub-grids have
    different indices than the ones displayed in this figure for the basic
    grid.

"""
function set_velocity_x_at_right_boundary_ghost_nodes!(
    vx::Float64,
    model::ModelData
)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bright = model.bcs.arrays.velocity.bright.array
    for i in 1:ynum+1
        bright[i,1] = vx
        bright[i,2] = 0.0
    end
    return nothing
end

end # module