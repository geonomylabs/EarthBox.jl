# Lava Flow Update Steps

[1] Calculate the reference eruption width ``W_{e}`` based on the characteristic magmatic crust 
height ``H_{mc(k)}`` for drainage basin ``k`` using equation [eq:eruption-width](@ref).

[2] Calculate minimum x-coordinate of the eruption domain ``x_{e,min}`` using equation 
[eq:x-eruption-min](@ref).

[3] Forecast the y-coordinate ``y_{e,fc}`` of the eruption by interpolating 
the y-coordinate of the topography grid ``(x_{topo}, y_{topo})`` at ``x_{e,fc}``, which is
calculated from [eq:x-eruption-fc](@ref) using the reference eruption width ``W_{e}``
and the minimum x-coordinate of the eruption domain ``x_{e,min}``.

[4] Calculate the characteristic volume per flow ``V_{flow}^{char}`` from equation
[eq:char-flow-volume](@ref) that uses the forecasted y-coordinate of the eruption ``y_{e,fc}``
to determine subaerial or submarine conditions.

[5] Calculate the number of flows per model time step ``N_{flow}`` from [eq:num-flows](@ref)
using the characteristic volume per flow ``V_{flow}^{char}`` and the total volume of volcanic
material extruded from the drainage basin ``V_{vol(k)}`` that was calculated during melt extraction.

[6] Calculate the volume of volcanic material per flow ``V_{flow}`` from equation 
[eq:flow-volume](@ref) using the number of flows per model time step ``N_{flow}`` and the
total volume of volcanic material extruded from the drainage basin ``V_{vol(k)}``.

[7] Define the initial y-coordinate of topography grid ``y_{topo, init}`` by setting it equal
to the current y-coordinate of the topography grid ``y_{topo}`` (this will be used when the compaction
section is added to the text).

[8] Initialize the integrated lava flow thickness ``H_{flow,tot}`` to zero.

[9] Run [Lava Flow Loop](@ref).

[10] Apply compaction correction to the y-coordinate of topography grid ``y_{topo}`` using 
an approach similar to that used to account for sediment compaction during the sediment
transport model but with the thickness of newly deposited lava flows ``H_{flow,tot}`` replacing
the sediment thickness. Pre-existing sedimentary and sticky markers are adjusted to account for 
compaction. The initial y-coordinate of the topography grid ``y_{topo,init}`` is used
to calculate the compaction correction.

## Lava Flow Loop

**For Each Lava Flow Do:**

[1] Calculate the x-coordinate of the eruption site ``x_{e}`` for the current flow as a 
normally distributed random variable using equations [eq:random-eruption-location](@ref) and 
[eq:eruption-location-constraints](@ref), the minimum x-coordinate of the eruption domain 
``x_{e,min}`` and the eruption width ``W_{e}``.

[2] Determine residual flow thickness ``H_{res}`` from equation [eq:residual_flow_thickness](@ref)
using the forecasted y-coordinate of the eruption ``y_{e,fc}`` and the y-coordinate of sea level 
``y_{sl}`` to determine subaerial or submarine conditions.

[3] Decimate the topography grid ``(x_{topo}, y_{topo})`` to obtain a grid with a lower 
resolution ``(x_{topo}', y_{topo}')`` that will be used for lava flow calculations using a 
cellular automata approach.

[4] Calculate a the decimated grid spacing ``\Delta x'``.

[5] Calculate the number of pulses per flow ``N_{pulse}`` from equation [eq:num_pulses](@ref)
using the flow volume ``V_{flow}``, the decimated grid spacing ``\Delta x'``, and the residual
flow thickness ``H_{res}``.

[6] Calculate pulse thickness, ``H_{pulse}``, from equation [eq:pulse_thickness](@ref) using 
the flow volume ``V_{flow}``, the number of pulses per flow ``N_{pulse}``, and the decimated grid 
spacing ``\Delta x'``.

[7] Initialize lava flow thickness ``H_{flow}'`` to zero on the decimated grid.

[8] Run [Lava Flow Pulse Loop](@ref).

[9] Calculate the updated thickness of the lava flow from multiple pulses on main grid 
``H_{flow}`` by interpolating ``H_{flow}'``.

[10] Update the y-coordinate of topography ``y_{topo}`` with flow thickness from equation 
[eq:topo-flow-update](@ref) using ``H_{flow}`` after copying ``y_{topo}`` to ``y_{topo}^o``. 

[11] Update total laval flow thickness integrating multiple flow events, ``H_{flow,tot}``, 
using the following equation:

###### eq:total-laval-flow-thickness
```math
    H_{flow,tot} = H_{flow,tot}^o + H_{flow}
```

where ``H_{flow,tot}^o`` is the previous total lava flow thickness.


## Lava Flow Pulse Loop

**For Each Lava Pulse ``p`` Do:**

[1] Calculate the index of the eruption node ``i_e`` on the decimated grid from 
[eq:eruption-index](@ref) using the x-coordinate of the eruption ``x_{e}`` and the 
decimated grid spacing ``\Delta x'``.

[2] Sort indices of the decimated grid ``(x_{topo}', y_{topo}')`` based on distance 
from the eruption node ``i_e``.

[3] Add the pulse thickness ``H_{pulse}`` to the flow thickness grid ``H_{flow}'`` at the 
eruption node ``i_e`` using equation [eq:initialize-pulse-thickness](@ref).

[4] Run [Cellular Automata Lava Flow Update Steps](@ref).


## Cellular Automata Lava Flow Update Steps

**While Maximum Thickness Difference ``>`` Tolerance Do for Each Node ``i``:**

- Update lava flow thickness ``H_{flow}'`` using cellular automata approach
whereby the flow thickness for each cell is updated starting from node ``i_e`` and 
radiating outward by allowing lava to flow from each cell with higher relative elevation to
adjacent cells with the lower relative elevation. 

!!! note "lava-flow-requirement"
    Lava is only allowed to flow to adjacent cells if thickness of lava is greater than the 
    residual thickness ``H_{res}``. 

!!! note "flow-thickness-update"
    Flow thickness is updated by iterating over grid nodes until the maximum thickness difference 
    between iterations falls below a user defined tolerance. The outflow thickness associated
    with flow from cell ``i`` to adjacent cells is calculated using equation 
    [eq:outflow-thickness](@ref) and the updated flow thickness is calculated using equation
    [eq:flow-thickness-update](@ref).