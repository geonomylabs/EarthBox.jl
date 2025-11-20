# Melt Extrusion and Lava Flow

Volcanic extrusion and lava flow for a given eruption time step is modeled using a cellular automata 
approach where the volume of magma available for extrusion in a melt drainage basin is divided into 
``N_{flow}`` flow events that extrude from a probabilistic eruption location ``x_{e}``. Flow events occur at 
an eruption interval ``\Delta t_e`` and involve the extrusion of all melt extracted during model time steps
within the eruption interval. At the beginning of each 
flow event, the current topography grid ``(x_{topo} \text{, } y_{topo})`` is decimated to a uniform 
grid ``(x_{topo}'\text{, } y_{topo}')`` with spacing ``\Delta x_{topo}'`` to improve the computational 
performance of the cellular automata calculation. Then a decimated flow thickness grid for the flow 
event, ``(x_{topo}'\text{, } H_{flow}')``, is initialized to zero thickness. Each flow is that divided 
into ``N_{pulse}`` pulses that laterally flow from the eruption location ``x_{e}`` to adjacent cells on 
the decimated topography grid. 

The volume of volcanic material per flow ``V_{flow}`` is calculated using the following equation:

###### eq:flow-volume
```math
V_{flow} = \frac{V_{vol,k}}{N_{flow}}
```

where ``N_{flow}`` is given by:

###### eq:num-flows
```math
N_{flow} =
\begin{cases}
    int\left( \frac{V_{vol,k}}{V_{flow}^{char}} \right) & \quad \text{for} \quad V_{vol,k} > V_{flow}^{char}\\
    1 & \quad \text{for} \quad \text{otherwise}
\end{cases}
```

where ``V_{vol,k}`` is the total volume of volcanic material available for extrusion from the drainage 
basin ``k`` calculated during the melt extraction steps within the eruption interval and ``V_{flow}^{char}`` 
is the characteristic volume of volcanic material per flow given by:

###### eq:char-flow-volume
```math
V_{flow}^{char} = 
\begin{cases}
    L_{sa} H_{res,sa}
        & \quad \text{for} \quad  y_{e,fc} \leq y_{sl}
            \quad \text{and} \\
    L_{sm} H_{res,sm}
        & \quad \text{for} \quad  y_{e,fc} > y_{sl}
\end{cases}
```

where ``L_{sa}`` is the characteristic subaerial flow length, ``H_{res,sa}`` is the 
subaerial residual flow thickness, ``L_{sm}`` is the characteristic submarine flow length, and
``H_{res,sm}`` is the submarine residual flow thickness, ``y_{e,fc}`` is the forecasted
y-coordinate of the eruption, and ``y_{sl}`` is the y-coordinate of sea level. The forecasted y-coordinate 
of the eruption ``y_{e,fc}`` in equation [Eq.](@ref eq:char-flow-volume) is determined by interpolating 
the y-coordinate from the decimated topography grid ``(x_{topo}' \text{, } y_{topo}')`` at the 
forecasted eruption location ``x_{e,fc}`` which is calculated using the following equation:

###### eq:x-eruption-fc
```math
x_{e,fc} = x_{e, min} + \frac{W_e}{2}
```

where ``x_{e, min}`` is the minimum x-coordinate of the eruption location give by:

###### eq:x-eruption-min
```math
x_{e, min} = x_{shallow,pm}' - \frac{W_e}{2}
```

where ``x_{shallow,pm}'`` is the average x-coordinate of the shallowest partially molten mantle 
marker within drainage basin ``k`` calculated during melt extraction, and ``W_e`` is the width 
of the eruption zone given by the following equation:

###### eq:eruption-width
```math
W_e = 
\begin{cases}
    W_{e,min}
        & \quad \text{for} \quad H_{mc} \leq H_{mc,min} \quad \text{and} \\
    W_{e,max}
        & \quad \text{for} \quad H_{mc} \geq H_{mc,max} \quad \text{and} \\
    W_{e,min} + 
        \frac{(H_{mc,k} - H_{mc,min})}{(H_{mc,max} - H_{mc,min})}(W_{e,max} - W_{e,min})
        & \quad \text{for} \quad H_{mc,min} < H_{mc} < H_{mc,max}
\end{cases}
```

where ``H_{mc,k}`` is the characteristic magmatic crust height for melt drainage basin ``k`` calculated 
during melt extraction, ``H_{mc,min}`` is the minimum characteristic magmatic crust height, ``H_{mc,max}`` 
is the maximum characteristic magmatic crust height, ``W_{e,min}`` is the minimum width of the 
eruption zone and ``W_{e,max}`` is the maximum width of the eruption zone. 

Each flow is erupted on the surface at eruption location ``x_{e}`` calculated as a normally distributed
random variable as described by the following equation:

###### eq:random-eruption-location
```math
x_e \sim \mathcal{N} \left( x_{e,min} + 0.5W_e, \, (0.25W_e)^2 \right)
```

subject to the constraint:

###### eq:eruption-location-constraints
```math
x_{e,min} \leq x_e \leq x_{e,min} + W_e \text{.}
```

If ``x_e`` falls outside this range, it is recalculated.

Each flow is erupted in a series of pulses where the number of pulses is:

###### eq:num_pulses
```math
N_{pulse} = \max\left[
    int\left(\frac{V_{flow}}{\left(\Delta x_{topo}' H_{res}\right)}\right)
    \text{, } 1
    \right]
```

where ``\Delta x_{topo}'`` is the grid spacing of a decimated topography grid used to perform
cellular automata calculations, and ``H_{res}`` is the residual flow thickness given by:

###### eq:residual_flow_thickness
```math
H_{res} = 
\begin{cases}
    H_{res,sa}
        & \quad \text{for} \quad  y_e \leq Y_{sl}
            \quad \text{and} \\
    H_{res,sm}
        & \quad \text{for} \quad  y_e > Y_{sl}
\end{cases}
```

The thickness of each pulse is calculated using the following equation:

###### eq:pulse_thickness
```math
H_{pulse} = \frac{V_{flow}}{N_{pulse}\Delta x_{topo}'}
```
A critical step in the cellular automata calculation is to sort the indices of the decimated topography
grid based on the distance from the eruption location ``x_{e}`` with grid index given by:

###### eq:eruption-index
```math
i_e = int\left(\frac{x_e}{\Delta x_{topo}'}\right) \text{.}
```

For a given pulse of volcanic material, the cellular automata calculation is initialized as:

###### eq:initialize-pulse-thickness
```math
H_{flow(i_e)}' = H_{flow(i_e)}'^{o} + H_{pulse} \text{.}
```

where ``H_{flow(i_e)}'^{o}`` is the thickness of the flow at index ``i_e`` prior to the current pulse.
For a given lava pulse, flow thickness ``H_{flow}'`` is updated for each cell by radiating output from 
index ``i_e``. Basaltic lava material flows from cells with higher relative elevation to adjacent cells 
with the lower relative elevation until the integrated thickness of the lava flow ``H_{flow}'`` equals 
the residual thickness ``H_{res}`` after which flow ceases to occur. During the iterative update procedure 
for a given pulse of the lava, the y-coordinates of topography, ``y_{topo(i)}''``, ``y_{topo(i-1)}''``, 
and ``y_{topo(i+1)}''`` are are related to current flow thickness ``H_{flow(i)}'`` as follows:

```math
\begin{split}
    y_{topo(i)}'' = y_{topo}(i)' - H_{flow(i)}' \\
    y_{topo(i-1)}'' = y_{topo}(i-1)' - H_{flow(i-1)}' \\
    y_{topo(i+1)}'' = y_{topo}(i+1)' - H_{flow(i+1)}'
\end{split}
```

and elevation differences that drive flow from cell ``i`` to adjacent cells are given by:

```math
\begin{split}
    \Delta E_{left(i)}' = \max(0 \text{, } y_{topo(i-1)}'' - y_{topo(i)}'') \\
    \Delta E_{right(i)}' = \max(0 \text{, } y_{topo(i+1)}'' - y_{topo(i)}'') \text{.}
\end{split}
```

where positive values of ``\Delta E_{left(i)}'`` and ``\Delta E_{right(i)}'`` indicate that flow will 
occur from cell ``i`` to adjacent cells ``i-1`` and ``i+1``, respectively, if current flow thickness
``H_{flow(i)}'`` is greater than the residual thickness ``H_{res}``. The total thickness of lava 
available to flow to adjacent cells is given by:

###### eq:total-outflow-thickness
```math
\Delta H_{out,total(i)}' = \frac{2}{3}
    \min\left[
            \max\left( 0 \text{, } H_{flow(i)}' - H_{res}\right)
            \text{, }
            \Delta E_{limit(i)}'
        \right]
```

where ``\Delta H_{out,total(i)}'`` is the total thickness of lava available to flow from cell ``i`` to
adjacent cells, and ``\Delta E_{limit(i)}'`` is given by:

###### eq:lava-delta-elevation-limit
```math
\Delta E_{limit(i)}' = 
\begin{cases}
    \left(\Delta E_{left(i)}' + \Delta E_{right(i)}'\right)/2
        & \quad \text{if} \quad
            \Delta E_{left(i)}' > 0 \quad \text{and} \quad 
            \Delta E_{right(i))}' > 0 \\
    \Delta E_{right(i)}'
        & \quad \text{else if} \quad \Delta E_{right(i)}' > 0 \\
    \Delta E_{left(i)}'
        & \quad \text{else} \quad \Delta E_{left(i)}' > 0 \text{.}
\end{cases}
```

Equation [Eq.](@ref eq:total-outflow-thickness) facilitates convergence of the cellular automata 
iterations by limiting the amount of inter-cellular flow and ensures that flow does not occur if
flow thickness is less than or equal to the residual flow thickness ``H_{res}``.

The total outflow thickness ``H_{out,total(i)}'`` is partitioned to adjacent cells using the following
equations:

###### eq:outflow-thickness
```math
\begin{split}
    \Delta H_{out,left(i)}' = min\left(
        \Delta H_{out,total(i)}'
        \frac{\Delta E_{left(i)}'}
            {\left( |\Delta E_{left(i)}'| + |\Delta E_{right(i)}'| \right)} \text{, } \Delta E_{left(i)}'
        \right) \\
    \Delta H_{out,right(i)}' = min\left(
        \Delta H_{out,total(i)}'
        \frac{\Delta E_{right(i)}'}
            {\left( |\Delta E_{left(i)}'| + |\Delta E_{right(i)}'| \right)} \text{, } \Delta E_{right(i)}'
        \right) \text{. }
\end{split}
```

where outflow thickness is limited by elevation differences between cells. The thickness of the flow
is updated for cell ``i`` and adjacent cells using the following equations:

###### eq:flow-thickness-update
```math
\begin{split}
    H_{flow(i-1)}' = H_{flow(i-1)}'^{o} + \Delta H_{out,left(i)}' \\
    H_{flow(i)}' = H_{flow(i)}'^{o} - \Delta H_{out,left(i)}' - \Delta H_{out,right(i)}' \\
    H_{flow(i+1)}' = H_{flow(i+1)}'^{o} + \Delta H_{out,right(i)}'
\end{split}
```

where ``H_{flow(i-1)}'^{o}``, ``H_{flow(i)}'^{o}``, and ``H_{flow(i+1)}'^{o}`` are the thickness of the flow
prior to the current cellular automata iteration. For each iteration of the cellular automata calculation,
equations [Eq.](@ref eq:flow-thickness-update) are applied to each cell in the decimated topography grid.
These iterations are repeated until the maximum thickness difference between iterations falls below a 
tolerance of ``10^{-4}``.

After the thickness of the flow on the decimated grid, ``H_{flow,i}'``, is updated for each pulse, flow
thickness is interpolated and added to the original un-decimated topography grid producing an 
updated flow thickness ``H_{flow,i}``. The topography grid at the end of a flow event is then updated as follows:

###### eq:topo-flow-update
```math
y_{topo,i} = y_{topo,i}^o - H_{flow,i}
```

where ``y_{topo,i}^o`` is y-coordinate of the topography from the prior flow event. Equation 
[Eq.](@ref eq:topo-flow-update) ensures that lava flow calculations are consistent with a dynamically 
evolving landscape impacted by the flow of volcanic material. 

The thickness of lava ``H_{flow,i}`` produced during a given flow event is added to ``H_{flow,tot,i}``, which is 
the integrated thickness of all flow events associated with a thermo-mechanical model time step.
The total flow thickness ``H_{flow,tot,i}`` is used to transform sticky markers to volcanic
material if ``y_m`` is greater than ``y_{topo,m}``, which is the interpolated y-coordinate of the
updated topography grid ``y_{topo,i}`` at x-coordinate ``x_m``. For cases where topography has been 
updated for both sedimentation and extrusive flow, sticky markers are transformed to volcanic 
sediment if ``y_m`` is greater than ``y_{topo,m} + H_{flow,tot,m}`` and to volcanic markers if 
``y_m`` is less than ``y_{topo,m} + H_{flow,tot,m}`` where ``H_{flow,tot,m}`` is the total thickness of lava 
flows interpolated at the x-coordinate of the marker ``x_m``. This approach assumes that sediment is
deposited prior to the eruption of volcanic material. The current volume of volcanic material available
for extrusion, ``V_{vol,k}``, is reset to zero at the end of each eruption interval.