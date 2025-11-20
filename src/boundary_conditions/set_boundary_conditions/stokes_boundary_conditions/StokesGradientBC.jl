module StokesGradientBC

import EarthBox.ModelDataContainer: ModelData

"""
    set_zero_vertical_vx_gradient_at_top_boundary!(model::ModelData)::Nothing

Define coefficients for staggered vx grid along top boundary.

# Updated arrays from group bcs.arrays.velocity
- btop.array::Matrix{Float64}

# Background
The x-component of staggered velocity vx is defined on a staggered grid with ghost 
nodes (Figure 1). For vx ghost nodes located at the top of the model and outside 
the main basic grid, the following equation is used to define vx:

vx[1,j] = btop[j,1] + vx[2,j]*btop[j,2]

where j is the x-index of the grid and the dimensions of vx are (ynum+1, xnum).

Free slip is implemented by setting vx ghost nodes vx[1,j] located outside the grid 
equal to the value of the next deepest staggered vx node vx[2,j] in the internal 
grid. For this case btop[j,1] = 0 and btop[j,2] = 1:

vx[1,j] = vx[2,j]

The array size of btop is xnum+1 since there are xnum+1 vy staggered nodes in the 
x-direction and the btop array is used to define coefficients for both vx and vy 
staggered grids. Since there are xnum vx staggered grid nodes in the x-direction 
there is and excess element in btop for the vx staggered grid.

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
function set_zero_vertical_vx_gradient_at_top_boundary!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    btop = model.bcs.arrays.velocity.btop.array
    for j in 1:xnum+1
        btop[j,1] = 0.0
        btop[j,2] = 1.0
    end
    return nothing
end

"""
    set_zero_vertical_vy_gradient_in_top_cells!(model::ModelData)::Nothing

Define coefficients staggered vy grid along top boundary.

# Updated arrays from group bcs.arrays.velocity
- btop.array::Matrix{Float64}

# Background
The y-component of velocity vy is defined on a staggered grid with ghost nodes 
(Figure 1). The following equation is used to define vy:

vy[1,j] = btop[j,3] + vy[2,j]*btop[j,4]

where j is the x-index of the grid and the dimensions of vy are (ynum, xnum+1).

Zero vertical gradient is implemented by setting vy nodes vy[1,j] along the top 
boundary boundary and outside the side boundaries of the basic grid equal to the 
value of the next deepest staggered vy node vy[2,j] in the internal grid. For this 
case btop[j,3] = 0 and btop[j,4] = 1:

vy[1,j] = vy[2,j]

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
function set_zero_vertical_vy_gradient_in_top_cells!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    btop = model.bcs.arrays.velocity.btop.array
    for j in 1:xnum+1
        btop[j,3] = 0.0
        btop[j,4] = 1.0
    end
    return nothing
end

"""
    set_zero_vertical_vx_gradient_at_bottom_boundary!(model::ModelData)::Nothing

Define coefficients for free slip along the bottom boundary for vx.

# Updated arrays from group bcs.arrays.velocity
- bbottom.array::Matrix{Float64}

# Background
The x-component of staggered velocity vx is defined on a staggered grid with ghost 
nodes (Figure 1). For vx ghost nodes located at the bottom of the model and outside 
the main basic grid, the following equation is used to define vx:

vx[ynum+1,j] = bbottom[j,1] + vx[ynum,j]*bbottom[j,2]

where j is the x-index of the grid and the dimensions of vx are (ynum+1, xnum).

Free slip is implemented by setting vx ghost nodes vx[ynum+1,j] located outside the 
grid equal to the value of the next shallowest staggered vx node vx[ynum,j] in the 
internal grid. For this case bbottom[j,1] = 0 and bbottom[j,2] = 1:

vx[ynum+1,j] = vx[ynum,j]

The array size of bbottom is xnum+1 since there are xnum+1 vy staggered nodes in 
the x-direction and the bbottom array is used to define coefficients for both vx 
and vy staggered grids. Since there are xnum vx staggered grid nodes in the 
x-direction there is and excess element in bbottom for the vx staggered grid.

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
function set_zero_vertical_vx_gradient_at_bottom_boundary!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    bbottom = model.bcs.arrays.velocity.bbottom.array
    for j in 1:xnum+1
        bbottom[j,1] = 0.0
        bbottom[j,2] = 1.0
    end
    return nothing
end

"""
    set_zero_vertical_vy_gradient_at_bottom_cells!(model::ModelData)::Nothing

Define coefficients for free slip along the bottom boundary for vy.

# Updated Arrays
- `model.bcs.arrays.velocity.bbottom.array`: Array((xnum+1,4), dtype=Float64)

# Background
The y-component of velocity vy is defined on a staggered grid with ghost nodes 
(Figure 1). The following equation is used to define vy:

vy[ynum+1,j] = bbottom[j,3] + vy[ynum,j]*bbottom[j,4]

where j is the x-index of the grid and the dimensions of vy are (ynum, xnum+1).

Free slip is implemented by setting vy nodes vy[ynum+1,j] along the bottom boundary 
boundary and outside the side boundaries of the basic grid equal to the value of 
the next shallowest staggered vy node vy[ynum,j] in the internal grid. For this 
case bbottom[j,3] = 0 and bbottom[j,4] = 1:

vy[ynum+1,j] = vy[ynum,j]

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
function set_zero_vertical_vy_gradient_at_bottom_cells!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    bbottom = model.bcs.arrays.velocity.bbottom.array
    for j in 1:xnum+1
        bbottom[j,3] = 0.0
        bbottom[j,4] = 1.0
    end
    return nothing
end

"""
    set_zero_lateral_vy_gradient_at_left_boundary!(model::ModelData)::Nothing

Define coefficients for free slip along the left boundary for vy.

# Updated arrays from group bcs.arrays.velocity
- bleft.array::Matrix{Float64}

# Background
The y-component of velocity vy is defined on a staggered grid with ghost nodes 
(Figure 1). The following equation is used to define vy:

vy[i,1] = bleft[i,3] + vy[i,2]*bleft[i,4]

where j is the x-index of the grid and the dimensions of vy are (ynum, xnum+1).

Free slip is implemented by setting vy at ghost nodes equal to the value of the 
next node toward the right in the internal grid.

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
function set_zero_lateral_vy_gradient_at_left_boundary!(model::ModelData)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bleft = model.bcs.arrays.velocity.bleft.array
    for i in 1:ynum+1
        bleft[i,3] = 0.0
        bleft[i,4] = 1.0
    end
    return nothing
end

"""
    set_zero_lateral_vy_gradient_at_right_boundary!(model::ModelData)::Nothing

Define coefficients for free slip along the right boundary for vy.

The y-component of velocity vy is defined on a staggered grid with ghost nodes 
(Figure 1). The following equation is used to define vy:

vy[i,xnum+2] = bright[i,3] + vx[i,xnum+1]*bright[i,4]

where j is the x-index of the grid and the dimensions of vy are (ynum, xnum+1).

Free slip is implemented by setting vy at ghost nodes equal to the value of the 
next node toward the left in the internal grid.

# Updated arrays from group bcs.arrays.velocity
- bright.array::Matrix{Float64}

# Description of the 2D Staggered Grid Along Right Boundary
                                j  4           5            6

                                    vx           vx           vx
                                i

                                1---+-----vy-----+-----vy-----+     vy
                                    |            |            |
                                    |            |            |
                                    vx    p'     vx    p'     vx
                                    |            |            |
                                    |            |            |
                                2---+-----vy-----+-----vy-----+     vy
                                    |            |            |
                                    |            |            |
                                    vx    p'     vx    p'     vx
                                    |            |            |
                                    |            |            |
                                3---+-----vy-----+-----vy-----+     vy
                                    |            |            |
                                    |            |            |
                                    vx     p'     vx    p'     vx
                                    |            |            |
                                    |            |            |
                                4---+-----vy-----+-----vy-----+     vy
                                    |            |            |
                                    |            |            |
                                    vx    p'      vx    p'     vx
                                    |            |            |
                                    |            |            |
                                5---+-----vy-----+-----vy-----+     vy

                                    vx           vx           vx

Figure 1: 2D staggered grid and numbering scheme used for discretizing Stokes and 
continuity equations. Indices i and j apply to basic grid nodes denoted with "+". 
Symbols vx and vy refer to staggered sub-grid nodes for the x- and y-components of 
velocity respectively. p' refers to the staggered sub-grid nodes for scaled pressure 
with nodes located in the center of basic grid cells. Note that staggered sub-grids 
have different indices than the ones displayed in this figure for the basic grid.
"""
function set_zero_lateral_vy_gradient_at_right_boundary!(model::ModelData)::Nothing
    ynum = model.grids.parameters.geometry.ynum.value
    bright = model.bcs.arrays.velocity.bright.array
    for i in 1:ynum+1
        bright[i,3] = 0.0
        bright[i,4] = 1.0
    end
    return nothing
end

end # module
