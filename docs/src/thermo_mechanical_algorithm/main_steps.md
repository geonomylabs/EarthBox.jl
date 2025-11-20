# Main Steps

The main steps of the thermomechanical algorithm are as follows:

1. Initialize grids and markers [Model Initialization Steps](@ref)
2. Run time loop [Time Loop Steps](@ref)


## Model Initialization Steps


### Staggered Grid Initialization Steps

[1] Initialize the x- and y-coordinates of the basic grid ``(x_{b,j}, y_{b,i})`` for
model dimensions ``(x_{size}, y_{size})`` where ``j`` ranges from 1 to 
``N_{x,b}`` and ``i`` ranges from 1 to ``N_{y,b}`` and ``N_{x,b}`` and 
``N_{y,b}`` are the number of grid points in the x- and y-directions, 
respectively. The convention used in this work is that y is the vertical 
direction with values increasing with depth and x is the horizontal direction.
Basic grid coordinates take into account local refinement.

[2] Initialize basic grid spacings ``\Delta x_{b,j}`` and ``\Delta y_{b,i}`` 
where ``j`` ranges from 1 to ``N_{x,b}-1`` and ``i`` ranges from 1 to ``N_{y,b}-1``.

[3] Initialize pressure grid x- and y-coordinates ``(x_{p,j}, y_{p,i})``
where ``j`` ranges from 0 to ``N_{x,b}-1`` and ``i`` ranges from 1 to ``N_{y,b}-1``.

[4] Initialize pressure grid spacings ``\Delta x_{p,j}`` and ``\Delta y_{p,i}``
where ``j`` ranges from 1 to ``N_{x,b}-2`` and ``i`` ranges from 1 to ``N_{y,b}-2``.

[5] Initialize y-direction coordinates of ``v_x``-grid  ``y_{v_x,i}`` where ``i``
ranges from 1 to ``N_{y,b}+1``.

[6] Initialize y-direction spacing of ``v_x``-grid ``\Delta y_{v_x,i}`` where ``i``
ranges from 1 to ``N_{y,b}``.

[7] Initialize x-direction coordinates of ``v_y``-grid ``x_{v_y,j}`` where ``j``
ranges from 1 to ``N_{x,b}+1``.

[8] Initialize x-direction spacing of ``v_y``-grid ``\Delta x_{v_y,j}`` where ``j``
ranges from 1 to ``N_{x,b}``.

### Initialize Topography Grid

The topography grid is initialized using the following steps:

[1] Initialize Eulerian topography grid coordinates ``(x_{t,j}, y_{t,j})`` 
where ``j`` ranges from 1 to ``N_{t}`` and ``N_{t}`` is the number of
topography grid points with uniform spacing ``\Delta x_{t}``.
[2] Initialize total volcanic extrusion thickness grid 
``(x_{t,j}, H_{t,j}^{vol})`` with ``H_{t,j}^{vol} = 0``.


### Marker Initialization Steps

[1] Initialize marker coordinates ``(x_m, y_m)`` using equations [eq:marker-coordinates](@ref).

[2] Initialize marker composition ``C_m`` and marker temperature ``T_m`` based on initial conditions.

[3] Initialize initial friction angle for marker ``\theta'_m{^o}`` using the initial friction
angle of the material associated with the marker composition ``C_m``.

[4] Initialize initial cohesion for marker ``\sigma_{c,m}^o`` using the initial cohesion
of the material associated with the marker composition ``C_m``.

[5] Initialize marker parameters:

```math
\begin{split}
    \text{Pressure: } P_m = 0 \\
    \text{Deviatoric normal stress: } \sigma_{(xx)m}' = 0 \\
    \text{Deviatoric shear stress: } \sigma_{(xy)m}' = 0 \\
    \text{Deviatoric normal strain rate: } \dot{\epsilon}_{xx(m)}' = 0 \\
    \text{Deviatoric shear strain rate: } \dot{\epsilon}_{xy(m)}' = 0 \\
    \text{Melt Fraction: } M_m = 0 \\
    \text{Extracted Melt Fraction: } M_{(extracted)m} = 0 \\
    \text{Extractable Melt Fraction} M_{(extractable)m} = 0\\
    \text{Strain: } \epsilon_m = 0 \\
    \text{Strain-rate Ratio: } R_{sr}(m) = 1.0
\end{split}
```

[6] Initialize marker density ``\rho_m`` based on ``C_m``, ``T_m`` and ``P_m``.

[7] Initialize marker thermal conductivity ``k_m`` based on ``C_m``, and ``T_m``.

[8] Initialize marker heat capacity ``C_{p(m)}`` based on ``C_m``, and ``T_m``.


## Time Loop Steps
[1] Set model time step $\Delta t$ equal to visco-elastic time step $\Delta t_{ve}$.

[2] Execute [Pre-solver Steps](@ref) including updating transient 
boundary conditions, updating melt fraction and melt-related marker composition, 
updating rock properties, updating frictional-plastic and creep rheologies, and 
interpolating marker parameters to staggered grid arrays.

[3] Create marker and grid output files at model time $t$.

[4] Execute [Stokes-solver Steps](@ref) including solving the 
visco-elasto-plastic Stokes-continuity equations on the staggered grid, adjusting time step 
$\Delta t$ using the updated velocity solution to ensure maximum cell displacement limits are 
not exceeded, updating marker stress taking grid and sub-grid changes and updating marker
strain rate, marker strain rate ratio and marker pressure using updated grid velocity and 
pressure solutions. 

[5] Execute [Heat-solver Steps](@ref) including updating marker 
temperature for grid and sub-grid thermal diffusion using the updated time step $\Delta t$.

[6] Execute [Advection-solver Steps](@ref) using adjusted time step 
$\Delta t$ including advecting markers and rotating stress tensors using a 4-th order Runge 
Kutta scheme, updating marker total and plastic strain and correcting temperature in the sticky 
air/water layer.

[7] Execute [Post-solver Steps](@ref) including updating melt extraction,
updating topography for deformation, sediment transport and lava flow, transforming
marker composition for surface processes, updating sea level using isostasy, updating the next 
eruption time, recycling markers that have exited the model domain, and updating marker 
serpentinization.
