# Staggered Grid Structure

An orthogonal 2D staggered grid composed of a basic grid and three staggered grids is used to 
discretize governing equations using a conservative finite difference formulation described by 
[gerya2010](@citet) ([Fig.](@ref fig:staggered-grid)). The basic grid, denoted by subscript
``b``, encompasses the main model domain ``(x_{size}, y_{size})`` with ``N_{x,b}`` nodes in the 
``x``-direction and ``N_{y,b}`` nodes in the ``y``-direction. The x- and y-coordinates of the basic 
grid are ``(x_{b,j}, y_{b,i})`` where ``j`` is in the ``x``-direction ranging from ``1`` to ``N_{x,b}`` 
and ``i`` is in the ``y``-direction ranging from ``1`` to ``N_{y,b}``. Both ``y_{b,i}`` and ``i`` 
increase with depth. The basic grid spacings in the ``x``- and ``y``-directions are denoted by 
``\Delta x_{b,j}`` where ``j`` ranges from ``1`` to ``N_{x,b}-1`` and ``\Delta y_{b,i}`` where ``i`` 
ranges from ``1`` to ``N_{y,b}-1``, respectively. Temperature unknowns ``T_{(i,j)_b}``, deviatoric 
shear stress ``\sigma'_{xy(i,j)_b}``, deviatoric shear strain rate ``\dot \epsilon'_{xy(i,j)_b}``, 
and visco-plastic shear viscosity ``\eta_{vp(i,j)_b}`` are defined on the nodes of the basic grid. 

## Grid Indexing

Basic grid cells and unknowns are indexed with respect to the upper-left node of 
basic-grid cells that have a global cell index defined as follows [Fig.](@ref fig:staggered-grid):

###### eq:basic-cell-index
```math
I_{cell}  = (j-1)(N_{y,b}-1) + i
```

where ``i`` and ``j`` are indices of the upper-left basic grid node associated with a cell. Unknowns 
``v_x``, ``v_y``, and ``P`` are numbered with the increasing cell index ``I_{cell}`` as follows:

###### eq:global-indices-of-stokes-unknowns
```math
\begin{split}
    I_{v_x} = 3I_{cell} - 2 & \quad \text{for x-velocity unknown} \\
    I_{v_y} = 3I_{cell} - 1 & \quad \text{for y-velocity unknown} \\
    I_{P} = 3I_{cell} & \quad \text{for pressure unknown}
\end{split}
```

Temperature unknowns ``T`` are numbered using the global index ``I_{T}``, which is defined as follows:

###### eq:global-indices-of-heat-unknowns
```math
I_{T} = (j-1)N_{y,b} + i
```

where ``i`` and ``j`` are indices of the basic grid [Fig.](@ref fig:staggered-grid). 

![Staggered grid](images/staggered_grid_4x4.png)
##### fig:staggered-grid
*Staggered grid used to discretize governing equations. Nodes of 
the basic grid are shown as black circles. Nodes of the velocity-x grid are denoted by
red symbols, nodes of the velocity-y grid are denoted by light blue symbols, and
nodes of the pressure grid are denoted by green symbols. Node index pairs for the 
basic, velocity-x, velocity-y and pressure grids are denoted with subscripts ``b``, ``v_x``, ``v_y`` 
and $p$, respectively. The global index for temperature unknowns ``I_T`` is denoted
with bold font, and global cell index ``I_{cell}`` used to discretize the Stokes-continuity 
equations are denoted by bold font in parentheses. The global index of Stoke-continuity unknowns
is shown with numbers inside of colored symbols. Colored circles represent normal unknowns where 
finite-difference stencils are applied. Ghost nodes located outside of the main 
model domain, which are used to approximate derivatives in the finite difference 
stencils, are shown as squares. Boundary nodes are shown with diamonds with thin 
outlines. Boundary unknown nodes are denoted with diamonds with thick outlines.
The unknowns defined at boundary unknown nodes have prescribed values leading to large matrix
coefficients equal to one and right-hand side terms equal prescribed velocity values.*