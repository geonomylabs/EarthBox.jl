# Algorithm for Solving the Heat Advection-Diffusion Equation

The Eulerian-Lagrangian marker-in-cell approach from [gerya2010](@citet) is used to solve equation 
[eq:energy](@ref) whereby markers are advected using a 4th-order Runge-Kutta scheme and marker 
temperature ``T_m`` is updated using the conductive temperature solution on Eulerian grid nodes and 
temperature changes associated with subgrid thermal diffusion:

**For each time step do:**

[1] Interpolate current marker temperature ``T_m`` to basic grid nodes ``T_{(i,j)_b}`` using 
equation [eq:bilinear_interp2grid](@ref).

[2] Calculate purely conductive temperature solution over time interval ``\Delta t`` on basic
grid ``T_{solu(i,j)_b}`` by solving equation [eq:energy](@ref) without advective terms and 
using current grid temperature ``T_{(i,j)_b}`` as the initial condition. An adaptive time-stepping 
scheme is applied whereby the time step is reduced if the change in temperature exceeds a threshold.

[3] Calculate the conductive temperature change on the Eulerian grid using the following equation:

###### eq:temperature_change
```math
\Delta T_{(i,j)_b} = T_{solu(i,j)_b} - T_{(i,j)_b}
```

[4] Update marker temperature ``T_m`` taking grid and sub-grid temperature changes into account 
[Subgrid Thermal Diffusion Steps](@ref).

[5] Advect markers using a 4th-order Runge-Kutta scheme.

[6] Advance model time ``t`` by ``\Delta t``.


## Subgrid Thermal Diffusion Steps

[1] Calculate the initial nodal temperature on markers ``T_{nodal,m}`` by interpolating 
``T_{(i,j)_b}`` back to markers using equation [eq:bilinear-interp2markers](@ref).

[2] Calculate the total sub-grid temperature difference between interpolated nodal temperatures
and current marker temperatures using the following equation:

###### eq:subgrid-temp-difference
```math
\Delta T_{sgt,m} = T_{nodal,m} - T_m
```

[3] Calculate the relaxed sub-grid temperature difference on markers using the following equation:

###### eq:subgrid-temp-relaxed
```math
\Delta T_{sg,m} = \Delta T_{sgt,m}
\left(1 - exp \left(-D_{sg}\frac{\Delta t}{t_{sg, m}}\right)\right)
```

where ``D_{sg}`` is the sub-grid diffusion coefficient, ``\Delta t`` is the time step and ``t_{sg}`` 
is the thermal diffusion time scale given by:

###### eq:thermal-diffusion-timescale
```math
t_{sg, m} = \frac{\rho_m C_{p,m}}{k_m\left(\frac{2}{\Delta x^2} + \frac{2}{\Delta y^2}\right)}
```

where ``\rho_m`` is the density of the marker, ``C_{p,m}`` is the specific heat capacity of the 
marker, ``k_m`` is the thermal conductivity of the marker, ``\Delta x`` and ``\Delta y`` are the 
grid spacings in the ``x`` and ``y`` directions, respectively.

[4] Calculate the relaxed sub-grid temperature difference on Eulerian grid 
``\Delta T_{sg(i,j)_b}`` by interpolating ``\Delta T_{sg,m}`` using equation 
[eq:bilinear_interp2grid](@ref).

[5] Calculate the remaining temperature change  ``\Delta T_{r(i,j)_b}`` on Eulerian nodes by 
removing the relaxed sub-grid temperature difference using the following equation:

###### eq:remaining-temp-change
```math
\Delta T_{r(i,j)_b} = \Delta T_{(i,j)_b} - \Delta T_{sg(i,j)}
```

[6] Interpolate ``\Delta T_{r(i,j)_b}`` to markers using equation 
[eq:bilinear-interp2markers](@ref) to obtain ``\Delta T_{r,m}``, which is the nodal 
temperature change with sub-grid effects removed.

[7] Update marker temperature ``T_m`` taking into account both sub-grid effects
and the remaining temperature change from the basic grid using the following equation:

###### eq:marker-temp-update
```math
T_m = T_m + \Delta T_{sg,m} + \Delta T_{r,m}
```

