module NoSlipBC

using EarthBox.ModelDataContainer: ModelData
using ..VelocityBC

function set_no_slip_bc_at_top_boundary!(model::ModelData)::Nothing
    set_vx_for_no_slip_at_top_boundary!(model)
    VelocityBC.set_velocity_y_at_top_boundary_ghost_nodes!(0.0, model)
    return nothing
end

function set_no_slip_bc_at_bottom_boundary!(model::ModelData)::Nothing
    set_vx_for_no_slip_at_bottom_boundary!(model)
    VelocityBC.set_velocity_y_at_bottom_boundary_ghost_nodes!(0.0, model)
    return nothing
end

function set_no_slip_bc_at_left_boundary!(model::ModelData)::Nothing
    set_vy_for_no_slip_at_left_boundary!(model)
    VelocityBC.set_velocity_x_at_left_boundary_ghost_nodes!(0.0, model)
    return nothing
end

function set_no_slip_bc_at_right_boundary!(model::ModelData)::Nothing
    set_vy_for_no_slip_at_right_boundary!(model)
    VelocityBC.set_velocity_x_at_right_boundary_ghost_nodes!(0.0, model)
    return nothing
end

function set_no_slip_with_shear_bc_at_right_boundary!(
    vy_shear::Float64,
    model::ModelData
)::Nothing
    set_vy_for_no_slip_with_shear_at_right_boundary!(vy_shear, model)
    VelocityBC.set_velocity_x_at_right_boundary_ghost_nodes!(0.0, model)
    return nothing
end

""" Define coefficients for no slip along the top boundary for vx.

# Updated Arrays
- `model.bcs.arrays.velocity.btop.array::Matrix{Float64}`: Array((xnum+1,4))

# Background
The x-component of velocity vx is defined on a staggered grid with ghost nodes
(Figure 1). The following equation is used to define vx on ghost nodes located at
the top of the model and outside the main grid:

vx[1,j] = btop[j,1] + vx[2,j]*btop[j,2]

where j is the x-index of the staggered grid and the dimensions of vx are
(ynum+1, xnum).

The no-slip condition is implemented by setting the x-velocity in the ghost nodes
equal to the inverse of the x-velocity in the next deepest vx node in the
internal grid:

vx[1,j] = -vx[2,j]

where btop[j,1] = 0.0 and btop[j,2] = -1.0. This condition ensures that the
x-velocity at the top boundary is zero.

# Description of the 2D Staggered Grid Along Top Boundary
```
   j   1            2            3             4           5            6
       vx           vx           vx           vx           vx           vx
i
1 vy     +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
         |            |            |            |            |            |
         |            |            |            |            |            |
         vx    p'     vx    p'     vx    p'     vx    p'     vx    p'     vx
         |            |            |            |            |            |
         |            |            |            |            |            |
2 vy      +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
         |            |            |            |            |            |
```

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and
continuity equations. Indices i and j apply to basic grid nodes denoted with "+".
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components
of velocity respectively. p' refers to the staggered sub-grid nodes for scaled
pressure with nodes located in the center of basic grid cells. Note that
staggered sub-grids have different indices than the ones displayed in this figure
for the basic grid.
"""
function set_vx_for_no_slip_at_top_boundary!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    btop = model.bcs.arrays.velocity.btop.array
    for j in 1:xnum+1
        btop[j,1] = 0.0
        btop[j,2] = -1.0
    end
    return nothing
end

""" Define coefficients for no slip along the bottom boundary for vx.

# Updated Arrays
- `model.bcs.arrays.velocity.bbottom.array::Matrix{Float64}`: Array((xnum+1,4))

# Background
The x-component of velocity vx is defined on a staggered grid with ghost nodes.
The following equation is used to define vx:

vx[ynum+1,j] = bbottom[j,1] + vx[ynum,j]*bbottom[j,2]

where j is the x-index of the staggered grid and the dimensions of vx are
(ynum+1, xnum).

The no-slip condition is implemented by setting the x-velocity in the ghost nodes
equal to the inverse of the x-velocity in the next shallowest vx node in the
internal grid:

vx[ynum+1,j] = -vx[ynum,j]

where bbottom[j,1] = 0.0 and bbottom[j,2] = -1.0. This condition ensures that
the x-velocity at the bottom boundary is zero.

# Description of the 2D Staggered Grid Along Bottom Boundary
```
   j   1            2            3             4           5            6
i
         |            |            |            |            |            |
3 vy      +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
         |            |            |            |            |            |
         |            |            |            |            |            |
         vx    p'     vx    p'     vx    p'     vx    p'     vx    p'     vx
         |            |            |            |            |            |
         |            |            |            |            |            |
4 vy      +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy
         |            |            |            |            |            |
         |            |            |            |            |            |
         vx    p'     vx    p'     vx    p'     vx    p'     vx    p'     vx
         |            |            |            |            |            |
         |            |            |            |            |            |
5 vy      +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy


         vx           vx           vx           vx           vx           vx
```

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and
continuity equations. Indices i and j apply to basic grid nodes denoted with "+".
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components
of velocity respectively. p' refers to the staggered sub-grid nodes for scaled
pressure with nodes located in the center of basic grid cells. Note that
staggered sub-grids have different indices than the ones displayed in this figure
for the basic grid.
"""
function set_vx_for_no_slip_at_bottom_boundary!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    bbottom = model.bcs.arrays.velocity.bbottom.array
    for j in 1:xnum+1
        bbottom[j,1] = 0.0
        bbottom[j,2] = -1.0
    end
    return nothing
end

""" Define coefficients for no slip along the left boundary for vy.

# Updated Arrays
- `model.bcs.arrays.velocity.bleft.array::Matrix{Float64}`: Array((ynum+1,4))

# Background
The y-component of velocity vy is defined on a staggered grid with ghost nodes
(Figure 1). The following equation is used to define vy on ghost nodes located at
the left of the model and outside the main grid:

vy[i,1] = bleft[i,3] + vy[i,2]*bleft[i,4]

where i is the y-index of the staggered grid and the dimensions of vy are
(ynum, xnum+1).

No slip is implemented by setting vy at ghost nodes equal to the opposite value
of the next node toward the right in the internal grid

vy[i,1] = -vy[i,2]

where bleft[i,3] = 0.0 and bleft[i,4] = -1.0. This condition ensures that the
y-velocity at the left boundary is zero.

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
2 vy      +-----vy-----+-----vy-----+--
         |            |            |
         |            |            |
         vx    p'     vx    p'     vx
         |            |            |
         |            |            |
3 vy      +-----vy-----+-----vy-----+--
         |            |            |
         |            |            |
         vx    p'     vx    p'     vx
         |            |            |
         |            |            |
4 vy      +-----vy-----+-----vy-----+---
         |            |8           |
         |            |            |
         vx    p'     vx    p'     vx
         |            |            |
         |            |            |
```

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and
continuity equations. Indices i and j apply to basic grid nodes denoted with "+".
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components
of velocity respectively. p' refers to the staggered sub-grid nodes for scaled
pressure with nodes located in the center of basic grid cells. Note that
staggered sub-grids have different indices than the ones displayed in this figure
for the basic grid.
"""
function set_vy_for_no_slip_at_left_boundary!(model::ModelData)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bleft = model.bcs.arrays.velocity.bleft.array
    for i in 1:ynum+1
        bleft[i,3] = 0.0
        bleft[i,4] = -1.0
    end
    return nothing
end

""" Define coefficients for no slip with shear along the right boundary.

# Updated Arrays
model.bcs.arrays.velocity.bright.array:: Matrix{Float64}: (ynum+1,4)

# Background
The y-component of velocity vy is defined on a staggered grid with ghost
nodes. The following equation is used to define vy on ghost nodes located
at the right of the model and outside the main grid:

vy[i,xnum+1] = bright[i,2] + vx[i,xnum]*bright[i,3]

where i is the y-index of the staggered grid and the dimensions of vy are
(ynum, xnum+1).

No slip is implemented by setting vy at ghost nodes equal to the opposite
value of the next node toward the left in the internal grid:

vy[i,xnum+1] = -vy[i,xnum]

where bright[i,2] = 0.0 and bright[i,3] = -1.0. This condition ensures that
the y-velocity at the right boundary is zero.

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
function set_vy_for_no_slip_at_right_boundary!(model::ModelData)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bright = model.bcs.arrays.velocity.bright.array
    for i in 1:ynum+1
        bright[i,3] = 0.0
        bright[i,4] = -1.0
    end
    return nothing
end

""" Define coefficients for no slip with shear along the right boundary.

# Updated Arrays
model.bcs.arrays.velocity.bright.array:: Matrix{Float64}: (ynum+1,4)

# Background
The y-component of velocity vy is defined on a staggered grid with ghost
nodes. The following equation is used to define vy on ghost nodes located
at the right of the model and outside the main grid:

vy[i,xnum+1] = bright[i,2] + vx[i,xnum]*bright[i,3]

where i is the y-index of the staggered grid and the dimensions of vy are
(ynum, xnum+1).

No slip with shear is implemented by setting vy at ghost nodes equal to
the opposite value of the next node toward the left in the internal grid
plus the shear velocity:

vy[i,xnum+1] = vy_shear + -vy[i,xnum]

where bright[i,2] = vy_shear and bright[i,3] = -1.0. This condition ensures
that the y-velocity at the right boundary is equal to the shear velocity.

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
function set_vy_for_no_slip_with_shear_at_right_boundary!(
    vy_shear::Float64,
    model::ModelData
)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bright = model.bcs.arrays.velocity.bright.array
    for i in 1:ynum+1
        bright[i,3] = vy_shear
        bright[i,4] = -1.0
    end
    return nothing
end

end # module