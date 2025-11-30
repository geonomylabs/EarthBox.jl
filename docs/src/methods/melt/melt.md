# Melt generation, Transport and Emplacement

The models presented in this work use a linearized equilibrium melt fraction model based on the 
following equation:

###### eq:melt-fraction
```math
M_m = \frac{T_m - T_{solidus,m}}{T_{liquidus,m} - T_{solidus,m}}
```

where ``m`` is the index of the marker, ``M_m`` is the melt fraction for the marker, ``T_m`` is the 
temperature of the marker, ``T_{solidus,m}`` is the solidus temperature of the marker, and 
``T_{liquidus,m}`` is the liquidus temperature of the marker [gerya13](@cite). 

The the solidus temperature of markers with ultramafic composition is defined for anhydrous 
peridotite by the following equation:

###### eq:solidus-mantle
```math
T_{solidus}^{mantle} =
\begin{cases}
    - 5.1(P_{GPa})^2 + 132.9P_{GPa} + 1085.7
        & \text{for } P_{GPa} \leq 8 \text{ GPa }\\
    T_{solidus,K03}^{8GPa}
    + \frac{\left( T_{solidus20} - T_{solidus,K3}^{8GPa}\right)}{12}
        \left(P_{GPa} - 8\right)
            & \text{for } P_{GPa} > 8 \text{ GPa}
\end{cases}
```

where ``P_{GPa}`` is the pressure in GPa, ``T_{solidusK3}^{8GPa}`` is the solidus temperature in ``K`` at 
8 GPa according to the model of [katz03](@citet), and ``T_{solidus20}`` is the solidus temperature in ``K``
at 20 GPa set equal to 2250 ``K``. Anhydrous peridotite is used as water is quickly removed from the 
system during the early stages of melting. The liquidus temperature for ultramafic mantle markers is
defined by the following equation for anhydrous peridotite:

###### eq:liquidus-mantle
```math
T_{liquidus}^{mantle} =
\begin{cases}
    - 2(P_{GPa})^2 + 45.0P_{GPa} + 1780.0
        & \text{for } P_{GPa} \leq 8 \text{ GPa }\\
    T_{liquidus,K3}^{8GPa}
    + \frac{\left( T_{liquidus20} - T_{liquidus,K03}^{8GPa}\right)}{12}
        \left(P_{GPa} - 8\right)
            & \text{for } P_{GPa} > 8 \text{ GPa}
\end{cases}
```

where ``T_{liquidusK3}^{8GPa}`` is the liquidus temperature in ``K`` at 8 GPa according to the model of
[katz03](@citet), and ``T_{liquidus20}`` is the liquidus temperature in ``K`` at 20 GPa set equal to 2600 ``K``.

This study employs two distinct models to determine the solidus and liquidus temperatures for gabbro: 
(1) a standard gabbro model representing the composition of primary melt products derived from the 
mantle, and (2) a model for gabbroic materials that have experienced fractional crystallization 
during emplacement. The solidus and liquidus temperatures for markers that have normal gabbro 
composition are defined by the following equations [gerya2010](@citep):

###### eq:solidus-liquidus-gabbro
```math
\begin{split}
    T_{solidus}^{gabbro} = 1327.0 + 91P_{GPa} \text{, }\\
    T_{liquidus}^{gabbro} = 1423 + 105P_{GPa}
\end{split}
```

The solidus and liquidus temperatures for fractionated gabbro are defined the following equations:

###### eq:solidus-liquidus-fractionated-gabbro
```math
\begin{split}
    T_{solidus}^{frac} = T_{solidus}^{gabbro} + \Delta T_{solidus}^{frac}\\
    T_{liquidus}^{frac} = T_{liquidus}^{gabbro} + \Delta T_{liquidus}^{frac}
\end{split}
```

where ``\Delta T_{solidus}^{frac}`` and ``\Delta T_{liquidus}^{frac}`` are the temperature differences 
between the solidus and liquidus temperatures for normal gabbro and fractionated gabbro, respectively. 
We use the following values for fractionated shifts in solidus and liquidus temperatures based on the
work of [maclennan04](@citet): ``\Delta T_{solidus}^{frac} = 100`` ``K`` and 
``\Delta T_{liquidus}^{frac} = 300`` ``K``. 

Melt is assumed to migrate instantaneously to the top of the partially molten mantle domain and 
then flow horizontally toward local high points [kneller25](@cite). At the local maximum, melt is 
divided into two fractions: one that is immediately emplaced at the base of the crust and another that 
is transported to the surface. This simplified model of melt transport and emplacement is implemented 
independently within each melt drainage basin bounded by local minima along the top of the partially 
molten domain.

The extractable melt fraction for a marker ``m`` is defined by the following 
equation:

###### eq:extractable-melt-fraction
```math
M_{extractable,m} = M_{m} - M_{extracted,m}^o
```

where ``M_{m}`` is the total melt fraction for the marker as described in [Eq.](@ref eq:melt-fraction) and 
``M_{extracted,m}^o`` is the melt fraction already extracted. The extracted melt fraction updated using:

###### eq:extracted-melt-fraction
```math
M_{extracted,m} = M_{extracted,m}^o + M_{extractable,m}
    \quad \text{if} \quad M_m > 0 \quad \text{and} \quad M_{extractable,m} > 0
```
Marker density ``\rho_m`` is impacted by melting by depleting the rock matrix and through the presence of
melt in the pore space. The impact of depletion is modeled using the approach of [theunissen22](@citet) 
as described in following equation:

###### eq:density-melt-depletion
```math
\rho_{m} = 
\begin{cases}
    \rho_m^o - 3.8 M_{extracted,m} 100
        & \text{if } M_{extracted,m} < 0.075 \\
    \rho_m^o - \left(24.84 + 0.488M_{extracted,m}100\right)
        & \text{if } M_{extracted,m} \geq 0.075 \\
\end{cases}
```

where ``\rho_m^o`` is the density of the marker prior to the update step. The effect of melt in pore 
space on density modeled using the following equation:

###### eq:density-melt-pore-space
```math
\begin{split}
    \rho_{melt} = \rho_{ref,melt} + f_{eos}P_{Pa} \text{, }\\
    \rho_m = \rho_m^o(1 - M_{extractable,m}) + \rho_{melt} M_{extractable,m}
\end{split}
```

where ``\rho_m^o`` is the density of the marker prior to this update step, ``\rho_{ref,melt}`` is the 
reference density of the melt set equal to 2750 ``\frac{kg}{m^3}``, ``P_{Pa,m}`` is the pressure in ``Pa`` 
and ``f_{eos}`` is the linearized derivative of melt density with respect to pressure equal to 
``8.97\cdot10^{-8} \frac{kg}{m^3Pa}`` derived from the Birch-Murnagham equation-of-state for anhydrous 
basaltic melt using ``K_{to} = 20.8`` and ``K_{tp} = 4.6``
[lesher15, stopler81](@citep).

The total number of magma markers in drainage basin ``k`` if all extractable melt 
were consolidated into a single body is calculated as:

###### eq:total-magma-markers
```math
N_{magma,tot,k}
    = f_{extraction}
        \left( \displaystyle \sum_{m=0}^{N_{marker}} M_{extractable,m} 
            + R_{melt,k}^o
        \right)
```

subject to the conditions:

```math
M_m > 0 \text{, } \quad M_{extractable,m} > 0 \quad \text{and} \quad
x_{divide,k} < x_m < x_{divide(k+1)}
```

where ``f_{extraction}`` is the extraction efficiency, ``R_{melt,k}^o`` is the residual 
melt from the previous time step, ``x_{divide,k}`` and ``x_{divide(k+1)}`` are the x-coordinates of the 
left and right drainage divides for drainage basin ``k``, respectively, ``x_m`` is the x-coordinate of 
the marker. For the experiments presented in this work ``f_{extraction} = 0.99`` to account for
a small amount of melt retained in the mantle during transport. 

The number of magma markers available for emplacement or extrusion in drainage basin ``k`` for a given 
time step is:

###### eq:available-magma-markers
```math
N_{magma,k} = int(N_{magma,tot,k})
```

and the residual component used in the next time step is updated as:

###### eq:residual-melt
```math
R_{melt,k} = N_{magma,tot,k} - N_{magma,k} \text{.}
```

The total extractable magma volume in the drainage basic ``k`` is then:

###### eq:total-extractable-magma-volume
```math
V_{magma,k} = N_{magma,k} \Delta x_m \Delta y_m
```

where ``\Delta x_m`` and ``\Delta y_m`` are the average widths of the markers in the x and y directions,
respectively.

The height of new column of crust formed if ``V_{magma,k}`` is emplaced as a dike within a spreading 
plate is:

###### eq:char-magmatic-crust-height
```math
H_{mc,k}
    = \frac{V_{magma,k}}{v_{ext}\Delta t}
```

where ``v_{ext}`` is the full extension velocity and ``\Delta t`` is the model time step. This 
characteristic magmatic crust height, ``H_{mc,k}``, is used to as a proxy for the efficiency of melt
transport to the surface.

The number of volcanic markers that are available for extrusion on the surface as lava flows above 
drainage basin ``k`` is given by:

###### eq:volcanic-markers
```math
N_{volcanic,k} = int(f_{v}N_{magma,k})
```

where ``f_{v}`` is the extrusion efficiency factor. We assume that the efficiency of melt transport to
the surface correlates with the characteristic magmatic crust height ``H_{mc,k}`` as described by the
following equation:

###### eq:extrusion_efficiency
```math
f_{v} = 
    \begin{cases}
        f_{v,min}
            & \quad \text{for} \quad H_{mc,k} \leq H_{mc,min}
                \quad \text{and} \\
        f_{v,max}
            & \quad \text{for} \quad H_{mc,k} \geq H_{mc,max}
                \quad \text{and} \\
        f_{v,min} + 
            \frac{(H_{mc,k} - H_{mc,min})}{\Delta H_{mc}}\Delta f_v
            & \quad \text{for} \quad H_{mc,min} < H_{mc,k} < H_{mc,max}
    \end{cases}
```

where ``\Delta f_v`` is equal to ``f_{v,max} - f_{v,min}`` and ``\Delta H_{mc}`` 
is equal to ``H_{mc,max} - H_{mc,min}``. For this work, ``f_{v,min} = 0.03``, ``H_{mc,min} = 6000 m`` and
``H_{mc,max} = 17500 m``. Alternative scenarios for parameter ``f_{v(max)}`` are explored ranging from
0.1 to 0.5. The total volume of volcanic markers ``V_{vol,k}`` extruded on the surface for a given
time step is:

###### eq:volcanic-extrusion-volume
```math
V_{vol,k} = V_{vol,k}^o + N_{volcanic,k} \Delta x_m \Delta y_m \text{.}
```

where ``V_{vol,k}^o`` is the volume of material available for extrusion from previous model times 
steps within the current eruption interval. The number of magma markers ``N_{magma,k}`` can be updated
to account for volcanic material as follows:

###### eq:magma-markers-final
```math
N_{magma,k} = N_{magma,k} - N_{volcanic,k} \text{.}
```

Therefore, the volume of magma markers from equation that are emplaced at the base of the crust is:

###### eq:magma-volume-final
```math
V_{magma,k} = N_{magma,k} \Delta x_m \Delta y_m \text{.}
```

Gabbroic magma is injected at the base of the Moho within an injection zone of width ``W_{inj,k}`` that
is centered on the local maximum of the partially molten domain of the mantle. The width of the magma 
injection zone ``W_{inj,k}`` is calculated using the following equation:

###### eq:injection-width
```math
W_{inj,k} = 
\begin{cases}
    \frac{V_{magma,k}}{H_{inj,limit}} &
        \quad \text{for} \quad H_{inj} > H_{inj,limit}\\
    W_{inj}^{char} &
        \quad \text{for} \quad H_{inj} \leq H_{inj,limit}
\end{cases}
```

where ``W_{inj}^{char}`` is the characteristic injection width, ``H_{inj,limit}`` is the maximum
injection height limit, and ``H_{inj}`` is the injection height defined by:

###### eq:injection-height
```math
H_{inj} = \frac{V_{magma,k}}{W_{inj}^{char}} \text{.}
```
The injection zone is divided into subdomains that are selected using a normal probability 
distribution centered on the local maximum of the partially molten domain of the mantle. The 
injection of gabbroic magma is implemented by replacing the shallowest marker below the Moho within 
a selected subdomain with a marker that has the composition of gabbroic magma and emplacement 
temperature of 1200``^{\circ}``C. At the beginning of each time step, gabbroic magma markers
that are within 2km of the Moho are converted into fractionated gabbroic magma markers. This 
conversion accounts for the fractional crystallization that occurs during the emplacement of magma 
in the crust, including processes such as dike and sill injection and the flow of crystal mush in 
the so called gabbro glacier.

Markers undergo solidification based on their updated melt fraction ``M_m``. The
solidification process occurs as follows: 

1. Solidification is modeled by transforming purely molten markers that have 
   cooled below the liquidus temperature to a solid state.

2. Solid markers with temperatures between the solidus and liquidus temperatures are 
   transformed to a partially molten state.

3. Partially molten markers that have cooled below the solidus are transformed to 
   a solid state by updating.

If the melt fraction ``M_m`` from [Eq.](@ref eq:melt-fraction) falls below the extracted melt 
fraction ``M_{extracted,m}`` from [Eq.](@ref eq:extractable-melt-fraction) (i.e. ``M_{extractable} < 0``) 
the marker is transformed into a refractory state.