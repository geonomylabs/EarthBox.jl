# Hydrothermal Update Steps

**For each marker ``m`` in drainage basin ``k`` do:**

1. Check if [Hydrothermal Cooling Criteria](@ref) are satisfied.
2. If criteria are met, update [Effective Thermal Conductivity](@ref).

## Hydrothermal Cooling Criteria

[1] Marker composition ``C_m`` is normal or fractionated gabbro, normal
or fractionated gabbroic emplaced magma, basalt, mantle or continental crust.

[2] There is a water column above the marker indicated by ``y_{ml} > y_{sl}``
where ``y_{ml}`` is the y-coordinate of the mudline at the x-coordinate ``x_m`` of the 
marker and ``y_{sl}`` is y-coordinate of sea level.

[3] The marker is below sea level as indicated by ``y_m > y_{sl}`` where ``y_m`` is
the y-coordinate of the marker.

[4] Sediment thickness at the x-coordinate of the marker ``d_{sed,x_m}`` is below a 
threshold of 2500 m.

## Effective Thermal Conductivity

[1] Calculate the x-limits of the molten zone ``x_{molten}^{min}`` and ``x_{molten}^{max}``
using equation [eq:x-molten-limits](@ref) with the x-coordinate of the mid-point of the 
gabbroic molten zone ``x_{molten,k}`` and the width of the gabbroic molten zone ``W_{molten,k}`` 
where ``k`` is the drainage basin index.

[2] Calculate the distance from the molten zone ``D_{molten}`` using equation 
[eq:distance-from-molten-zone](@ref) with ``x_{molten}^{min}``, ``x_{molten}^{max}`` and 
the x-coordinate of the marker ``x_m``.

[3] Calculate the hydrothermal distance factor ``f_{hydro}`` using equation 
[eq:hydrothermal-distance-factor](@ref) with ``D_{molten}``.

[4] Calculate the distance-based Nusselt number ``Nu_{dist}`` using equation 
[eq:distance-based-nusselt](@ref) with ``f_{hydro}``.

[5] Calculate the effective Nusselt number ``Nu_{eff}`` using equation 
[eq:effective-nusselt](@ref) with ``Nu_{dist}``, the plastic strain rate 
``\dot \epsilon_{plastic,m}`` and plastic strain ``\epsilon_{plastic(m)}`` of the marker ``m``.

[6] Update marker thermal conductivity ``k_m`` using equation [eq:hydrothermal-k](@ref) with
the current thermal conductivity ``k_m`` as the unmodified thermal conductivity of the marker,
marker x-coordinate ``x_m``, marker y-coordinate ``y_m``, marker temperature ``T_m``, the submud 
depth of the marker ``d_{sm,m}``, sediment thickness at x-coordinate of marker ``d_{sed,x_m}``, 
y-coordinate of mudline at the x-coordinate of marker ``y_{ml}`` and the effective Nusselt 
number ``Nu_{eff}``.