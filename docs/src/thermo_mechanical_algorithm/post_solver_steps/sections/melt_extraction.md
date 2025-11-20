# Melt Extraction Update Steps

[1] Calculate melt drainage divides ([Drainage Steps](@ref))

[2] Update extractable melt fraction ``M_{extractable,m}`` for each marker ``m`` using equation 
[eq:extractable-melt-fraction](@ref) with melt fraction ``M_m`` and and extracted melt fraction
``M_{extracted,m}``.

[3] Update extracted melt fraction ``M_{extracted,m}`` for each marker using equation 
[eq:extracted-melt-fraction](@ref) after copying current extracted melt fraction to ``M_{extracted}^o``.

[4] Calculate melt volumetrics for each drainage basin ``k`` including number 
of magma markers ``N_{magma(k)}``, volume of extrusive volcanic markers ``V_{volcanic(k)}``,
and the width of the injection domain ``W_{inj(k)}`` ([Melt Volumetrics Steps](@ref)).

[5] Inject magma at the base of the Moho above melt focusing point in each drainage basin ``k`` 
([Magma Injection Steps](@ref)).

[6] Calculate the dimensions of molten domain for each drainage basin ``k`` including
the width of the molten zone ``W_{molten,k}``, the mid point of the molten zone
``x_{molten,k}``, and the depth of the molten zone at the mid point ``y_{molten,k}``.
The molten zone includes all markers within the drainage basin ``k`` that have a
composition of normal gabbroic magma, partially molten gabbro, fractionated gabbroic
magma and partially molten fractionated gabbro.

[7] Update marker composition ``C_m`` by transforming all mantle markers to a
refractory state if ``M_{extractable} < 0``. 


## Drainage Steps

[1] Calculate the coordinates ``(X_{topo}, Y_{pm})`` of the top of the 
partially molten zone in the mantle by searching for the shallowest 
partially molten mantle marker at each topography grid point.

[2] Apply a low-pass filter similar to [eq:low-pass-filter](@ref) to the 
topography grid ``Y_{topo}`` to remove high-frequency noise.

[3] Calculate x-coordinates, ``X_{divides}``, of ``N_{divides}`` drainage divides, 
defined as local low points along the top of the partially molten zone 
``(X_{topo}, Y_{pm}^{lowpass})``. For each grid point calculate the slope to 
the left ``s_{left}`` and the slope to the right ``s_{right}``. The right slope 
is defined as the next non-zero slope to the right of the grid point. If the 
left slope is positive and the right slope is negative, the midpoint between
the current grid point and the next grid point where the right slope is 
negative is a drainage divide. Note that the convention used here involves 
the y-coordinates increasing downward. Drainage divides are assumed at the 
left and right boundaries of the model. 


## Melt Volumetrics Steps

[1] Calculate the number of partially molten mantle markers.

[2] Calculate the total number of magma markers ``N_{magma,tot,k}`` in 
drainage basin ``k`` if all extractable melt were consolidated into a single 
body using equation [eq:total-magma-markers](@ref) with ``M_{extractable,m}``
and ``R_{melt,k}^o``.

[3] Calculate the number of magma markers available for emplacement or extrusion
``N_{magma,k}`` using equation [eq:available-magma-markers](@ref) with ``N_{magma,tot,k}``.

[4] Update the residual fractional melt ``R_{melt(k)}`` using equation
[eq:residual-melt](@ref) with ``N_{magma,tot,k}`` and ``N_{magma,k}``.

[5] Calculate the total extractable magma volume ``V_{magma,k}`` using equation
[eq:total-extractable-magma-volume](@ref) with ``N_{magma,k}``.

[6] Calculate the characteristic magmatic crust height ``H_{mc,k}`` using equation
[eq:char-magmatic-crust-height](@ref) with ``V_{magma,k}``.

[7] Calculate the number of volcanic markers that are available for extrusion
``N_{volcanic,k}`` using equation [eq:volcanic-markers](@ref) with ``N_{magma,k}``
and ``H_{mc,k}``.

[8] Calculate the volcanic extrusion volume ``V_{vol,k}`` for the current 
drainage basin ``k`` using equation [eq:volcanic-extrusion-volume](@ref) with
``N_{volcanic,k}``.

[9] Update the number of magma markers ``N_{magma,k}`` to account for
the number of marker that will be extruded ``N_{volcanic,k}`` using equation 
[eq:magma-markers-final](@ref).

[10] Calculate the volume of magma available for emplacement ``V_{magma,k}`` using
equation [eq:magma-volume-final](@ref) with ``N_{magma,k}``.

[11] Calculate the injection height ``H_{inj}`` using equation [eq:injection-height](@ref)
with ``V_{magma,k}``.

[12] Calculate the injection width ``W_{inj,k}`` for magma emplacement at
the base of the Moho taking into account the injection height limit
``H_{inj,limit}`` using equation [eq:injection-width](@ref) with ``V_{magma,k}`` and ``H_{inj}``. 


## Magma Injection Steps

**For ``N_{magma(k)}`` markers do:**
    
[1] Determine the coordinates ``(x_{shallow}, y_{shallow})`` of the
shallowest partially molten mantle marker in drainage basin ``k``.

[2] Determine the coordinates ``(x_{shallow}', y_{shallow}')`` of
the shallowest marker with index ``m_{shallow}`` below the Moho within
a subdomain of the injection domain centered on 
``(x_{shallow}, y_{shallow})`` with width ``W_{inj,k}``. The injection domain
is divided into ``N_{sub}`` subdomains that are selected based on a normal 
probability distribution.

[3] Update marker composition ``C_m`` by transforming marker with 
index ``m_{shallow}`` to a gabbroic magma marker. Adjust marker 
parameters for transformation to magma and injection conditions:

###### eq:magma-parameter-adjustment
```math
\begin{split}
    \sigma_{xx,m}' = \sigma_{xy,m}' =  \epsilon_m = 
    \epsilon_{plastic,m} = \dot{\epsilon}_{xx,m} = \dot{\epsilon}_{xy,m} = \\
    M_{m} = M_{extractable,m} = M_{extracted,m} = R_{serp} = 0\\
    \theta_m = \theta_{mat}\text{,} \quad R_{sr} = 1\text{,}
    \quad T_m = T_{injection}\\
\end{split}
```

where ``\theta_{mat}`` is the initial friction angle for the material 
associated with the marker and ``T_{injection}`` is the melt injection 
temperature.