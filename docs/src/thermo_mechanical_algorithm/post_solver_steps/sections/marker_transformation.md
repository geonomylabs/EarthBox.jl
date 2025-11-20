# Marker Transformation Steps

**For each marker ``m`` do:**

[1] Determine the marker y-coordinate of topography ``y_{topo,m}`` by interpolating the topography 
grid ``(x_{topo,i}, y_{topo,i})`` at x-coordinate of marker ``x_m``.

[2] Determine the marker extrusion thickness ``H_{flow,tot,m}`` by interpolating
the topography grid ``(x_{topo,i}, H_{flow,tot,i})`` at x-coordinate of the marker ``x_m``.

[3] Define the y-coordinate of the top of the newly deposited sedimentary layer ``y_{sed,m}``
assuming new volcanic flows were emplaced after sedimentary layers were deposited using the 
following equation:

###### eq:top-new_sediment
```math
y_{sed,m} = y_{topo,m} + H_{flow,tot,m}
```

[4] Transform markers for water-air interface:
1. If sticky-air markers are below sea level ``y_{sl}``, transform sticky-air markers to
    water markers.
2. If sticky-water markers are above sea level ``y_{sl}``, transform sticky-water markers
    to air markers.

[5] If marker composition is sticky-air or sticky-water and ``y_m > y_{topo,m}`` transform markers
to sedimentary markers if ``y_m > y_{sed,m}`` or transform sticky-air markers to volcanic markers if
``y_m \leq y_{sed,m}``.

[6] Apply erosion if marker composition is lithological (i.e. not sticky-air or sticky-water)
and ``y_m < y_{interp,m}``:
1. If ``y_m < y_{sl}`` transform markers to sticky-air markers.
2. If ``y_m \geq y_{sl}`` transform markers to sticky-water markers.