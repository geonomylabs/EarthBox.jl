# Post-solver Steps

[1] Calculate melt extraction ([Melt Extraction Update Steps](@ref)).

[2] Update topography grid ``(x_{topo,i}, y_{topo,i})`` for material displacement using updated 
velocity solution ([Topography Grid Update Steps](@ref)).

[3] Apply the sediment transport model with sediment compaction to update topography grid 
``(x_{topo,i}, y_{topo,i})`` ([Sediment Transport Update Steps](@ref)).

[4] Copy topography grid ``(x_{topo}, y_{topo})`` to ``(x_{topo}', y_{topo}')``, which
will be used to track the evolving topography associated with lava flows.

[5] Apply lava flow model ([Lava Flow Update Steps](@ref)).

[6] Update sea level ``y_{sl}`` assuming local isostatic equilibrium between reference
lithosphere column and the average pressure at the base of the model.

[7] Update next eruption time ``t_{eruption}`` to ``t_{eruption} + \Delta t_e`` if the model 
time ``t`` is greater than the current eruption time ``t_{eruption}``.

[8] Recycle markers that have exited the model domain by re-injecting rock markers at the bottom
of the model domain as asthenosphere, and by re-injecting sticky air-water markers at the top
of the model domain.

[9] Apply marker transformation for erosion, sedimentation and extrusion by updating marker
composition ``C_m`` ([Marker Transformation Steps](@ref)).

[10] Calculate the serpentinization ratio ``R_{serp,m}`` using equation 
[eq:serpentinization-ratio](@ref) with marker temperature ``T_m`` and plastic strain rate 
``\dot \epsilon_{plastic,m}`` after copying the current serpentinization ratio to ``R_{serp,m}^o``.
