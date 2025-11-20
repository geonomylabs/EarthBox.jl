# Melt Damage Update Steps

For each drainage basin ``k`` do

1. Calculate the central damage zone damage probability ``p_{damage}^{central}``using equation 
   [eq:central-damage-probability](@ref) with characteristic magmatic crust height 
   ``H_{mc,k}`` calculated during melt extraction.
2. Calculate the x-coordinate limits of the melt damage zone ``x_{damage}^{min}`` and 
   ``x_{damage}^{max}`` and the y-coordinate limit ``y_{damage}^{max}`` using equations
   [eq:melt-damage-zone-limits](@ref) with the average x-coordinate of the shallowest partially
   molten mantle marker ``x'_{shallow,pm,k}`` and average y-coordinate of the shallowest partially
   molten mantle marker ``y'_{shallow,pm,k}`` calculated during melt extraction.
3. Calculate the melt damage factor for each marker ([Melt Damage Marker Update](@ref)).

## Melt Damage Marker Update

For each marker ``m`` do:

1. Calculate the damage probability of the marker ``p_{damage,m}`` using equation 
   [eq:melt-damage-probability](@ref) with the x-coordinate of the marker ``x_m``,
   ``x'_{shallow,pm,k}``, ``p_{damage}^{central}``, ``x_{damage}^{min}`` and
   ``x_{damage}^{max}``.
2. Calculate a random number ``r`` between 0 and 1.
3. Assign the marker damage factor ``f_{damage,m}`` using equation 
   [eq:melt-damage-factor](@ref) with ``p_{damage,m}`` and ``r``.
