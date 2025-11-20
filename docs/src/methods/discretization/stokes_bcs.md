 # Stokes-continuity Boundary Conditions

Earthbox uses the approach of [gerya2010](@citet) for specifying boundary conditions  
for the Stokes-continuity equations on a Staggered grid whereby simple coefficients are used
that control type of boundary condition. 

## Stokes Boundary Conditions: Top and Bottom Boundaries

Boundary conditions for the x-component of velocity are implemented along the top and bottom 
boundaries using the following equations:

###### eq:vx-bc-top-bottom
```math
\begin{split}
    v_{x(1,j)_{v_x}} & = v_{x}^{top} + f_{v_x}^{top}v_{x(2,j)_{v_x}} \\
    v_{x(N_{y,b}+1,j)_{v_x}} & = v_{x}^{bot} + f_{v_x}^{bot}v_{x(N_{y,b},j)_{v_x}} \\
\end{split}
```

where ``j`` ranges from ``1`` to ``N_{x,b}``, ``f_{v_x}^{top}`` and ``f_{v_x}^{bot}`` are the 
boundary condition coefficients for the top boundary and bottom boundaries, respectively, and 
``v_{x}^{top}`` and ``v_{x}^{bot}`` are the prescribed velocity components at ghost nodes along 
the top and bottom boundaries, respectively. The coefficients ``f_{v_x}^{top}`` and ``f_{v_x}^{bot}`` 
control how the x-component of velocity at ghost nodes is related to the internal velocity and can 
have a value of -1, 0 or 1 depending on the type of boundary condition.

Boundary conditions for the y-component of velocity are implemented along the top and bottom 
boundaries using the following equations:

###### eq:vy-bc-top-bottom
```math
\begin{split}
    v_{y(1,j)_{v_y}} & = v_{y}^{top} + f_{v_y}^{top}v_{y(2,j)_{v_y}} \\
    v_{y(N_{y,b},j)_{v_y}} & = v_{y}^{bot} + f_{v_y}^{bot}v_{y(N_{y,b}-1,j)_{v_y}} \\
\end{split}
```

where ``j`` ranges from ``1`` to ``N_{x,b}+1``, ``f_{v_y}^{top}`` and ``f_{v_y}^{bot}`` are the 
boundary condition coefficients for the top boundary and bottom boundaries, respectively, and
``v_{y}^{top}`` and ``v_{y}^{bot}`` are the prescribed velocity components at boundary unknown nodes 
along the top and bottom boundaries, respectively. The coefficients ``f_{v_y}^{top}`` and 
``f_{v_y}^{bot}`` control how the y-component of velocity at boundary unknown nodes is related to 
the internal velocity and can have a value of -1, 0 or 1 depending on the type of boundary condition.

## Stokes Boundary Conditions: Left and Right Boundaries

Boundary conditions for the x-component of velocity are implemented along the left and right 
boundaries using the following equations:

###### eq:vx-bc-left-right
```math
\begin{split}
    v_{x(i,1)_{v_x}} & = v_{x}^{left} + f_{v_x}^{left}v_{x(i,2)_{v_x}} \\
    v_{x(i,N_{x,b})_{v_x}} & = v_{x}^{right} + f_{v_x}^{right}v_{x(i,N_{x,b}-1)_{v_x}} \\
\end{split}
```

where ``i`` ranges from ``1`` to ``N_{y,b}+1``, ``f_{v_x}^{left}`` and ``f_{v_x}^{right}`` are 
the boundary condition coefficients for the left boundary and right boundaries, respectively, 
``v_{x}^{left}`` and ``v_{x}^{right}`` are the prescribed velocity components at boundary unknown 
nodes along the left and right boundaries, respectively. The coefficients ``f_{v_x}^{left}`` and 
``f_{v_x}^{right}`` control how the x-component of velocity at boundary unknown nodes is related to 
the internal velocity and can have a value of -1, 0 or 1 depending on the type of boundary condition.

Boundary conditions for the y-component of velocity are implemented along the left and right 
boundaries using the following equations:

###### eq:vy-bc-left-right
```math
\begin{split}
    v_{y(i,1)_{v_y}} & = v_{y}^{left} + f_{v_y}^{left}v_{y(i,2)_{v_y}} \\
    v_{y(i,N_{x,b}+1)_{v_y}} & = v_{y}^{right} + f_{v_y}^{right}v_{y(i,N_{x,b})_{v_y}} \\
\end{split}
```

where ``i`` ranges from ``1`` to ``N_{y,b}``, ``f_{v_y}^{left}`` and ``f_{v_y}^{right}`` are 
the boundary condition coefficients for the left boundary and right boundaries, respectively, and 
``v_{y}^{left}`` and ``v_{y}^{right}`` are the prescribed velocity components at boundary unknown 
nodes along the left and right boundaries. The coefficients ``f_{v_y}^{left}`` and ``f_{v_y}^{right}`` 
control how the y-component of velocity at boundary unknown nodes is related to the internal velocity 
and can have a value of -1, 0 or 1 depending on the type of boundary condition.

## No-slip Boundary Conditions

For no-slip boundary conditions along all boundaries, boundary coefficients ``f_{v_x}^{top}``, 
``f_{v_x}^{bot}``, ``f_{v_y}^{left}`` and ``f_{v_y}^{right}`` are set equal to -1, and boundary
coefficients ``f_{v_y}^{top}``, ``f_{v_y}^{bot}``, ``f_{v_x}^{left}`` and ``f_{v_x}^{right}`` 
are set equal to 0. Additionally, prescribed velocities ``v_{x}^{top}``, ``v_{x}^{bot}``, 
``v_{y}^{top}``, ``v_{y}^{bot}``, ``v_{x}^{left}``, ``v_{x}^{right}``, ``v_{y}^{left}`` and 
``v_{y}^{right}`` are set to zero. These conditions ensure that velocity at ghost nodes has the 
opposite sign of the internal velocity along the boundary, making the velocity at the model 
boundary zero as described in the following equations:

###### eq:no-slip-bc
```math
\begin{split}
    v_{x(1,j)_{v_x}} & = -v_{x(2,j)_{v_x}} \\
    v_{y(1,j)_{v_y}} & = 0 \\
    v_{x(N_{y,b}+1,j)_{v_x}} & = -v_{x(N_{y,b},j)_{v_x}} \\
    v_{y(N_{y,b},j)_{v_y}} & = 0 \\
    v_{x(i,1)_{v_x}} & = 0 \\
    v_{y(i,1)_{v_y}} & = -v_{y(i,2)_{v_y}} \\
    v_{x(i,N_{x,b})_{v_x}} & = 0 \\
    v_{y(i,N_{x,b}+1)_{v_y}} & = -v_{y(i,N_{x,b})_{v_y}} \\
\end{split}
```

## Free-slip Boundary Conditions with Inflow/Outflow

For free slip with inflow/outflow boundary conditions, boundary coefficients ``f_{v_x}^{top}``,
``f_{v_x}^{bot}``, ``f_{v_y}^{left}`` and ``f_{v_y}^{right}`` are set equal to 1, and boundary
coefficients ``f_{v_y}^{top}``, ``f_{v_y}^{bot}``, ``f_{v_x}^{left}`` and ``f_{v_x}^{right}`` 
are set equal to 0. Additionally, prescribed velocities ``v_{x}^{top}``, ``v_{x}^{bot}``, 
``v_{y}^{left}`` and ``v_{y}^{right}`` are set to zero, and prescribed velocities ``v_{y}^{top}``, 
``v_{y}^{bot}``, ``v_{x}^{left}`` and ``v_{x}^{right}`` are set to the inflow/outflow velocity.
These conditions ensure that velocity parallel to the boundary at ghost nodes is equal to the 
internal velocity, and velocity perpendicular to the boundary at boundary nodes is equal 
to the inflow/outflow velocity as described in the following equations:

###### eq:free-slip-bc-inflow-outflow
```math
\begin{split}
    v_{x(1,j)_{v_x}} & = v_{x(2,j)_{v_x}} \\
    v_{y(1,j)_{v_y}} & = v_{y}^{top} \\
    v_{x(N_{y,b}+1,j)_{v_x}} & = v_{x(N_{y,b},j)_{v_x}} \\
    v_{y(N_{y,b},j)_{v_y}} & = v_{y}^{bot} \\
    v_{x(i,1)_{v_x}} & = v_{x}^{left} \\
    v_{y(i,1)_{v_y}} & = v_{y(i,2)_{v_y}} \\
    v_{x(i,N_{x,b})_{v_x}} & = v_{x}^{right} \\
    v_{y(i,N_{x,b}+1)_{v_y}} & = v_{y(i,N_{x,b})_{v_y}} \\
\end{split}
```

!!! tip "pure-free-slip"
    To implement standard free-slip condition, the inflow-outflow velocities are set to zero.

## Boundary Conditions for Extensional Models with Fixed Boundaries

Extensional models with fixed boundaries use free-slip boundary conditions with inflow/outflow along 
all boundaries as described in [Eq.](@ref eq:free-slip-bc-inflow-outflow). Prescribed outflow 
velocity along the side boundaries is set to ``0.5v_{ext}`` where ``v_{ext}`` is the full 
extension velocity. Prescribed inflow velocity along the top and bottom boundaries are defined 
so that the net inflow/outflow volume is zero over the entire model domain and the net 
inflow/outflow volume is zero within the sticky layer. The prescribed inflow velocity along 
the top boundary ``v_y^{top}`` is calculated as:

###### eq:inflow-velocity-top
```math
v_y^{top} = \frac{v_{ext}\left(H_{sticky,l} + H_{sticky,r}\right)}{x_{size}}
```

where ``H_{sticky,l}`` and ``H_{sticky,r}`` are the thickness of the sticky layer on the left 
and right sides of the model domain, and ``x_{size}`` is the width of the model domain. The 
prescribed inflow velocity along the bottom boundary ``v_y^{bot}`` is calculated as:

###### eq:inflow-velocity-bottom
```math
v_{y}^{bot} = v_y^{top} -2v_{ext}\frac{y_{size}}{x_{size}}
```

where ``y_{size}`` is the height of the model domain. Inflow velocities ``v_{y}^{top}`` and 
``v_{y}^{bot}`` are recalculated at each time step to account for changes in the sticky layer 
thickness due to crustal thickness changes and isostatic adjustment.