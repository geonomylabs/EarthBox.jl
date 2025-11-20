# Algorithm for Solving the Visco-elasto-plastic Stokes-Continuity Equations

The frictional-plastic failure model used in this work introduces a strong non-linearity where the
visco-plastic viscosity ``{\eta_{vp,m}}`` is dependent on deviatoric stress and 
deviatoric stress is in turn dependent on the visco-plastic viscosity as described in 
equations [eq:visco-plastic-viscosity](@ref) and [eq:second-invariant-elastic-stress](@ref). 
EarthBox uses a node-based Picard iteration method to iteratively solve the non-linear Stokes-continuity 
equations on the staggered grid. 

The deviatoric stress tensor is updated on markers that are advected with the velocity field. Updating 
the deviatoric stress tensor on markers requires accounting for three components of stress change:

1. deviatoric grid stress changes from equation [Eq.](@ref eq:visco-elasto-plastic-stress-changes)
   obtained on the Eulerian staggered grid by solving the non-linear Stokes-continuity equation 
   with visco-elasto-plastic stress forecast
2. deviatoric stress changes occurring at a subgrid scale on individual markers
3. deviatoric stress changes associated with the rotation of markers during advection

The following steps are used to solve equations [eq:x-stokes](@ref), [eq:y-stokes](@ref), 
and [eq:continuity](@ref) using the visco-elasto-plastic constitutive law described in equations 
[eq:xx-stress](@ref) and [eq:xy-stress](@ref), and update the deviatoric stress tensor on 
markers for the three components of stress change described above.

**For each time step ``t`` do:**

[1] Set model time step ``\Delta t`` equal to visco-elastic time step ``\Delta t_{ve}``.

[2] Update marker failure properties including cohesion ``\sigma_{c,m}`` and friction angle ``\theta_m``
including the effects of randomized initial friction angle ``\theta_m^o``, strain 
weakening and melt damage.

[3] Update non-linear creep viscosity ``\eta_{creep,m}`` including the effects of temperature, pressure,
stress.

[4] Initialize visco-plastic viscosity ``\eta_{vp,m}`` with ``\eta_{creep,m}`` and update for the effects
of partial melt.

[5] Interpolate Stokes-continuity parameters from markers to staggered grids using equation 
[eq:bilinear_interp2grid](@ref):
```math
\begin{split}
    \text{Density (inclusive): } \rho_m \rightarrow \rho_{(i,j)_b} \\
    \text{Cohesion (exclusive): } C_m \rightarrow C_{(i,j)_b} \\
    \text{Friction angle (exclusive): } \theta_m \rightarrow \theta_{(i,j)_b} \\
    \text{Creep viscosity (exclusive): } \eta_{creep,m} \rightarrow \eta_{creep{(i,j)_b}} \\
    \text{Visco-plastic Shear Viscosity (exclusive): } \eta_{vp,m} \rightarrow \eta_{vp(i,j)_b} \\
    \text{Visco-plastic Normal Viscosity (exclusive): } \eta_{vp,m} \rightarrow \eta_{vp(i,j)_p} \\
    \text{Shear modulus on basic nodes (exclusive): } \mu_m \rightarrow \mu_{(i,j)_b} \\
    \text{Shear modulus on pressure nodes (inclusive): } \mu_m \rightarrow \mu_{(i,j)_p} \\
    \text{Deviatoric shear stress (exclusive): } \sigma_{xy,m}' \rightarrow \sigma_{xy(i,j)_{b1}}' \\
    \text{Deviatoric normal stress (inclusive): } \sigma_{xx,m}' \rightarrow \sigma_{xx(i,j)_{p1}}' \\
\end{split}
```

[6] Calculate updated solutions for x-component of velocity ``v_{x(i,j)_{vx}}``, y-component of velocity 
``v_{y(i,j)_{vy}}`` and pressure ``P_{(i,j)_p}`` using a node-based Picard iteration loop 
([Picard Loop](@ref)) with iterative updates to visco-plastic shear viscosity on the basic grid
``\eta_{vp(i,j)_b}`` that take into account plastic yielding and visco-plastic normal viscosity on pressure 
grid ``\eta_{vp(i,j)_p}`` calculated using a harmonic average of updated shear visco-plastic viscosity 
on basic nodes.

[7] Update model time step ``\Delta t`` ensuring that displacement is less than a user specified
fraction of the staggered grid spacing based on current velocity solutions ``v_{x(i,j)vx}`` and 
``v_{y(i,j)vy}`` and grid spacings ``\Delta x_{(i,j)vx}`` and ``\Delta y_{(i,j)vy}``.

[8] Calculate marker visco-plastic viscosity ``\eta_{vp,m}`` by interpolating visco-plastic viscosity 
from basic grid ``\eta_{vp(i,j)_b}`` and taking plastic yielding into account
([Marker Visco-plastic Viscosity](@ref)).

[9] Forecast deviatoric stresses ``\sigma_{xx(i,j)_{p2}}'`` and ``\sigma_{xy(i,j)_{b2}}'`` using 
equation [Eq.](@ref eq:deviatoric-stress-forecast) with updated deviatoric strain rates 
``\dot{\epsilon}_{xx(i,j)_p}'`` and ``\dot{\epsilon}_{xy(i,j)_b}'``, updated visco-plastic 
viscosity ``\eta_{vp(i,j)_p}`` and ``\eta_{vp(i,j)_b}`` and updated model time step ``\Delta t``.

[10] Update grid deviatoric stress changes using the following equations:

###### eq:visco-elasto-plastic-stress-changes
```math
\begin{split}
\Delta \sigma_{xx(i,j)_p}' = \sigma_{xx(i,j)_{p2}}' - \sigma_{xx(i,j)_{p1}}' \\
\Delta \sigma_{xy(i,j)_b}' = \sigma_{xy(i,j)_{b2}}' - \sigma_{xy(i,j)_{b1}}'
\end{split}
```

[11] Update marker deviatoric stress components ``\sigma_{xx,m}'`` and ``\sigma_{xx,m}'`` by applying
subgrid-stress diffusion steps ([Subgrid-stress Diffusion Steps](@ref)).

[12] Advect markers and calculate marker angular velocity ``\omega_m`` using the 4-th order 
Runge-Kutta method and angular velocity on the basic grid ``\omega_{(i,j)_b}`` given by:

###### eq:angular-velocity
```math
\omega_{(i,j)_b} =
    0.5\left(
        \frac{(v_{y(i,j+1)vy} - v_{y(i,j)vy})}{ \Delta x_{vy(j)} }
        - \frac{(v_{x(i+1,j)vx} - v_{x(i,j)vx})}{\Delta y_{vx(i)}}
        \right)
```

[13] Rotate deviatoric stress components on markers using the following equations:

###### eq:stress-rotation
```math
\begin{split}
\sigma_{xx,m}' = 
    \sigma_{xx,m}'^o\left(\cos(\omega_m \Delta t)^{2} - \sin(\omega_m \Delta t)^{2}\right)
    - \sigma_{xy,m}'^o \sin(2\omega_m \Delta t
    ) \text{, }\\
\sigma_{xy,m}' = 
    \sigma_{xx,m}'^o \sin(2\omega_m \Delta t) + \sigma_{xy,m}'^o \cos(2\omega_m \Delta t)
\end{split}
```

where ``\sigma_{xx,m}'^o`` and ``\sigma_{xy,m}'^o`` are the deviatoric marker stress components 
before advection.

[14] Advance model time ``t`` by ``\Delta t``.


## Picard Loop

**While Picard criterion ``R_{L2}`` ``>`` ``Tolerance`` and ``N_{iterations} < N_{max}`` do:**

[1] Calculate the visco-plastic viscosity on pressure grid ``\eta_{vp(i,j)_p}`` as a harmonic
    average of the visco-plastic viscosity ``\eta_{vp(i,j)_b}`` on surrounding basic nodes:

###### eq:harmonic-viscosity-average
```math
\eta_{vp(i,j)_p} = \frac{4}{
    1/\eta_{vp(i,j)_b} + 1/\eta_{vp(i-1,j)_b} 
    + 1/\eta_{vp(i,j-1)_b} + 1/\eta_{vp(i-1,j-1)_b}
    }
```

[2] Obtain new solution vector ``S^{new}`` and updated solutions for x-component of velocity 
    ``v_{x(i,j)_{vx}}``, y-component of velocity ``v_{y(i,j)_{vy}}`` and pressure ``P_{(i,j)_p}`` 
    by solving equations [eq:x-stokes](@ref), [eq:y-stokes](@ref), and [eq:continuity](@ref) 
    using current ``\eta_{vp(i,j)_b}`` and ``\eta_{vp(i,j)_p}``.

[3] Update the convergence ``R_{L2}`` criterion defined as:

###### eq:picard_criterion
```math
R_{L2} = \frac{
    \sqrt{\displaystyle\sum_{i} (S_{v_{xy},i}^{new} - S_{v_{xy},i}^{old})^2}
    }{
    \sqrt{ \displaystyle\sum_{i} (S_{v_{xy},i}^{new})^2}
    }
```

where ``S_{v_{xy}}^{old}`` is the old velocity solution vector and ``S_{v_{xy}}^{new}`` is the 
new velocity solution vector.

[4] Calculate deviatoric strain rate components ``\dot{\epsilon}_{xx(i,j)_p}'`` and 
``\dot{\epsilon}_{xy(i,j)_b}'`` using updated velocity solutions ``v_{x(i,j)vx}`` and 
``v_{y(i,j)vy}`` and the following equation:

###### eq:deviatoric-strain-rate
```math
\begin{split}
\dot{\epsilon'}_{xx(i,j)_p} = \frac{1}{2} \left( 
    \frac{
        v_{x(i+1,j+1)vx} - v_{x(i+1,j)vx}
        }{\Delta x_{j_b}} 
    - \frac{
        v_{y(i+1,j+1)vy} - v_{y(i,j+1)vy}
        }{\Delta y_{i_b}} 
\right) \quad \text{and} \\
\dot{\epsilon'}_{xy(i,j)_b} = \frac{1}{2} \left( 
    \frac{
        v_{x(i+1,j)vx} - v_{x(i,j)vx}
        }{\Delta y_{i_p}} 
    + \frac{
        v_{y(i,j+1)vx} - v_{y(i,j)vx}
        }{\Delta x_{j_p}}
\right) \text{, }
\end{split}
```

[5] Forecast deviatoric stresses ``\sigma_{xx(i,j)_{p2}}'`` and ``\sigma_{xy(i,j)_{b2}}'`` using 
updated deviatoric strain rates ``\dot{\epsilon}_{xx(i,j)_p}'`` and ``\dot{\epsilon}_{xy(i,j)_b}'``, 
updated visco-plastic viscosity ``\eta_{vp(i,j)_p}`` and ``\eta_{vp(i,j)_b}`` with the 
following equations:

###### eq:deviatoric-stress-forecast
```math
\begin{split}
    \sigma_{xx(i,j)_{p2}}' = 
        2 \eta_{vp(i,j)_p} \dot{\epsilon}_{xx(i,j)_p}'
            \frac{ \mu_{(i,j)_b} \Delta t }{ \mu_{(i,j)_b} \Delta t + \eta_{vp(i,j)_p} } 
        + \sigma_{xx(i,j)_{p1}}' 
            \frac{ \eta_{vp(i,j)_p} }{ \mu_{(i,j)_b} \Delta t + \eta_{vp(i,j)_p} } 
                \quad \text{and} \\
    \sigma_{xy(i,j)_{b2}}' =
        2 \eta_{vp(i,j)_b} \dot{\epsilon}_{xy(i,j)_b}'
            \frac{ \mu_{(i,j)_b} \Delta t }{ \mu_{(i,j)_b} \Delta t + \eta_{vp(i,j)_b} }
        + \sigma_{xy(i,j)_{b1}}'
            \frac{ \eta_{vp(i,j)_b} }{ \mu_{(i,j)_b} \Delta t + \eta_{vp(i,j)_b} }
\end{split}
```

[6] Update visco-plastic viscosity on basic grid ``\eta_{vp(i,j)_b}`` for plastic yielding,
and update yield state ``\chi_{(i,j)_b}`` 
([Visco-plastic Viscosity Update for Yielding on Basic Grid](@ref)).


## Visco-plastic Viscosity Update for Yielding on Basic Grid

**For Each Basic Node ``(i,j)`` Do:**

[1] Calculate the second invariant of deviatoric stress ``\sigma_{{II}(i,j)_b}'`` using updated 
deviatoric stresses ``\sigma_{xx(i,j)_{b2}}'`` and ``\sigma_{xy(i,j)_{b2}}'`` as described by the 
following equation:

###### eq:second-invariant-stress_basic_nodes
```math
\begin{split}
\sigma_{xx,avg(i,j)_{b2}}' =
    \frac{1}{4}
    \left(
        \sigma_{xx(i,j)_{p2}}' + \sigma_{xx(i+1,j)_{p2}}'
        + \sigma_{xx(i,j+1)_{p2}}' + \sigma_{xx(i+1,j+1)_{p2}}'
    \right) \text{,} \\
\sigma_{II(i,j)_b}' = \sqrt{
    \left(
        \left( \sigma_{xy(i,j)_{b2}}' \right)^2
        + \left( \sigma_{xx,avg(i,j)_{b2}}'\right)^2
    \right)
    }
\end{split}
```

[2] Compute the second invariant of deviatoric stress for purely elastic stress buildup
``\sigma_{II, elastic(i,j)_b}'`` using the following equation:

###### eq:elastic_stress
```math
\sigma_{II, elastic(i,j)_b}' = \frac{\mu_{(i,j)_b} \Delta t + \eta_{vp(i,j)_b}}{\eta_{vp(i,j)_b}}
    \sigma_{II(i,j)_b}'
```

where ``\mu_{(i,j)_b}`` is the shear modulus at basic nodes and ``\Delta t`` is the model time step.

[3] Interpolate pressure from staggered pressure grid ``P_{(i,j)_p}`` to basic pressure grid 
``P_{(i,j)_b}``:

###### eq:pressure_basic
```math
P_{(i,j)_b} = \frac{1}{4}\left(P_{(i,j)_p} + P_{(i+1,j)p} + P_{(i,j+1)p} + P_{(i+1,j+1)p}\right)
```

[4] Compute the Drucker-Prager yield stress on basic grid ``\sigma_{\text{yield}(i,j)_b}`` using
the following equation:

###### eq:yield_stress_basic_nodes
```math
\sigma_{yield(i,j)} =
\begin{cases}
    \sigma_{c(i,j)_b} \cos(\theta_{(i,j)_b}) + \sin(\theta_{(i,j)_b})\left(P_{(i,j)_b} - P_f\right) 
        & \quad \text{for} P_{(i,j)_b} \geq P_f\\
    \sigma_{c(i,j)_b} \cos(\theta_{(i,j)_b}) & \quad \text{for} P_{(i,j)_b} < P_f
\end{cases}
```

where ``\sigma_{yield(i,j)}`` is the second invariant of the deviatoric stress tensor at basic 
grid node at yield, ``\sigma_{c(i,j)_b}`` is the cohesion, ``\theta_{(i,j)_b}`` is the friction 
angle and ``P_f`` is the fluid pressure.
            
[5] Update visco-plastic viscosity ``\eta_{vp(i,j)_b}`` for plastic failure using the following
equation:

###### eq:visco_elasto_plastic_viscosity_update
```math
\begin{split}
\eta_{vp(i,j)_b} = 
    \begin{cases}
    \eta_{creep(i,j)_b} 
            & \text{for} \quad \sigma_{II,elastic(i,j)_b}' \leq \sigma_{yield(i,j)_b} \quad \\
    \min \left( \eta_{yield(i,j)_b}
            \text{, } \eta_{creep(i,j)_b} \right)
            & \text{for} \quad \sigma_{II,elastic(i,j)_b}' > \sigma_{yield(i,j)_b}
    \end{cases}
\end{split}
```

where ``\eta_{creep(i,j)_b}`` is the composite creep viscosity interpolated from markers and 
``\eta_{yield(i,j)_b}`` is the visco-plastic viscosity at the yield state as given by:

###### eq:yield_viscosity
```math
\eta_{yield(i,j)_b} = \mu_{(i,j)_b} \Delta t
        \frac{\sigma_{yield(i,j)_b}}{\sigma_{II,elastic(i,j)_b}' - \sigma_{yield(i,j)_b}}
```

[6] Update the plastic yielding state on basic nodes ``\chi_{(i,j)_b}`` using the following 
    equation:

###### eq:plastic_yielding
```math
\chi_{(i,j)_b} = 
\begin{cases}
    1 & \text{if } \sigma_{II,elastic(i,j)_b}' > \sigma_{yield(i,j)_b} 
        \text{ and } \eta_{yield(i,j)_b} < \eta_{creep(i,j)_b} \\
    0 & \text{otherwise}
\end{cases}
```

where a value of 1 indicates that the basic node has undergone plastic failure and a value of 0
indicates that the basic node has not undergone plastic failure.


## Marker Visco-plastic Viscosity

**For each marker ``m`` do:**

[1] Determine if plastic yielding has occurred for any basic node associated with the cell that
    contains marker ``m`` using the plastic failure criterion ``\chi_{(i,j)_b}`` defined on the 
    basic grid.
    
[2] Initialize the visco-plastic viscosity of the marker ``\eta_{vp,m}`` to the creep viscosity
    ``\eta_{creep(m)}``.

[3] If plastic yielding has occurred on a node associated with marker ``m`` calculate the 
    visco-plastic viscosity of the marker ``\eta_{vp,m}`` using a harmonic average limited by the 
    flow viscosity ``\eta_{creep,m}`` as follows:

###### eq:harmonic_viscosity_average_yield
```math
\begin{split}
\eta_{yield,m} =
    \left(
        \chi_{(UL)_b}w_{UL,m} + \chi_{(LL)_b}w_{LL,m}
        + \chi_{(UR)_b}w_{UR,m} + \chi_{(LR)_b}w_{LR,m}
    \right)\\
    \left(
    \frac{\chi_{(UL)_b}w_{UL,m}}{\eta_{vp(UL)_b}} +
    \frac{\chi_{(LL)_b}w_{LL,m}}{\eta_{vp(LL)_b}} +
    \frac{\chi_{(UR)_b}w_{UR,m}}{\eta_{vp(UR)_b}} +
    \frac{\chi_{(LR)_b}w_{LR,m}}{\eta_{vp(LR)_b}}
    \right)^{-1} \text{, } \\
\eta_{vp,m} = \min\left(\eta_{yield,m} \text{, } \eta_{creep,m}\right)
\end{split}
```

where ``UL = (i_{ul,m}, j_{ul,m})``, ``LL = (i_{ul,m}+1, j_{ul,m})``, ``UR = (i_{ul,m}, j_{ul,m}+1)``,
and ``LR = (i_{ul,m}+1, j_{ul,m}+1)`` are the basic nodes indices associated with the cell 
containing marker ``m`` and ``w_{UL,m}``, ``w_{LL,m}``, ``w_{UR,m}``, and ``w_{LR,m}`` are the 
associated node-marker weights calculated using equation eq:node_marker_weights
and ``i_{ul,m}`` and ``j_{ul,m}`` are the indices of the upper-left basic grid node of the cell
containing marker ``m``.


## Subgrid-stress Diffusion Steps

[1] Calculate initial nodal stresses ``\sigma_{xx,nodal,m}'`` and ``\sigma_{xy,nodal,m}'`` for each
    marker ``m`` by interpolating ``\sigma_{xx(i,j)_{p1}}'`` and ``\sigma_{xy(i,j)_{b1}}'`` to 
    markers using equation [eq:bilinear-interp2markers](@ref).

[2] Calculated the total subgrid stress differences ``\Delta \sigma_{xx,sgt,m}'`` and 
    ``\Delta \sigma_{xy,sgt,m}'`` for each marker ``m`` using the following equations:

###### eq:subgrid_stress_difference
```math
\begin{split}
\Delta \sigma_{xx,sgt,m}' & = \sigma_{xx,nodal,m}' - \sigma_{xx,m}' \\
\Delta \sigma_{xy,sgt,m}' & = \sigma_{xy,nodal,m}' - \sigma_{xy,m}' \text{.}
\end{split}
```

[3] Calculate the relaxed subgrid stress differences ``\Delta \sigma_{xx,sg,m}'`` and 
    ``\Delta \sigma_{xy,sg,m}'`` using the following equations:

###### eq:subgrid_stress_relaxed
```math
\begin{split}
\Delta \sigma_{xx,sg,m}' = \Delta \sigma_{xx,sgt,m}'
    \left(1 - \exp \left({-D_{sg} \frac{\Delta t}{t_{max}}}\right)\right) \\
\Delta \sigma_{xy,sg,m}' = \Delta \sigma_{xy,sgt,m}'
    \left(1 - \exp \left(-D_{sg} \frac{\Delta t}{t_{max}}\right)\right)
\end{split}
```

where ``D_{sg}`` is the sub-grid diffusion coefficient and ``t_{max}`` is the Maxwell time
defined as:

###### eq:maxwell-time
```math
t_{max} = \frac{\eta_{vp,m}}{\mu_m} \text{.}
```

[4] Calculate the relaxed sub-grid stress differences on the staggered grids 
    ``\Delta \sigma_{xx,sg(i,j)_p}'`` and ``\Delta \sigma_{xy,sg(i,j)_b}'`` by interpolating the 
    relaxed subgrid stress differences from markers ``\Delta \sigma_{xx,sg,m}'`` and 
    ``\Delta \sigma_{xy,sg,m}'`` to the staggered grid using the bilinear interpolation method 
    described in equation [eq:bilinear_interp2grid](@ref).

[5] Calculate the remaining stress changes on staggered grid nodes by removing the interpolated
    relaxed sub-grid stress differences using the following equations:

###### eq:remaining_stress_change_on_grid
```math
\begin{split}
\Delta \sigma_{xx,vepr(i,j)_p}' = 
    \Delta \sigma_{xx(i,j)_p}' - \Delta \sigma_{xx,sg(i,j)_p}' \text{, }\\
\Delta \sigma_{xy,vepr(i,j)_b}' = 
    \Delta \sigma_{xy(i,j)_b}' - \Delta \sigma_{xy,sg(i,j)_b}'
\end{split}
```

[6] Calculate the remaining visco-elasto plastic stress changes
``\Delta \sigma_{xx,vepr,m}'`` and ``\Delta \sigma_{xy,vepr,m}'`` for each marker ``m`` by 
interpolating the remaining stress changes on staggered grid nodes 
``\Delta \sigma_{xx,vepr(i,j)_b}'`` and ``\Delta \sigma_{xy,vepr(i,j)_b}'`` to markers using 
equation [Eq.](@ref eq:bilinear-interp2markers).

[7] Update marker deviatoric stress components ``\sigma_{xx,m}'`` and ``\sigma_{xx,m}'`` using
    the following equations:

###### eq:subgrid_stress_update
```math
\begin{split}
\sigma_{xx,m}' = 
    \sigma_{xx,m}'^{o} + \Delta \sigma_{xx,sg,m}' 
    + \Delta \sigma_{xx,vepr,m,}' \text{,}\\
\sigma_{xy,m}' = 
    \sigma_{xy,m}'^{o} + \Delta \sigma_{xy,sg,m}' 
    + \Delta \sigma_{xy,vepr,m}'
\end{split}
```

where ``\sigma_{xx,m}'^{o}`` and ``\sigma_{xy,m}'^{o}`` are the previous marker deviatoric stress
components.