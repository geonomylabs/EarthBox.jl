# Hydrothermal Circulation

Hydrothermal circulation reduces temperatures in the lithosphere at a regional scale in areas where 
locally hot rock associated magmatic intrusions interacts with seawater. Hydrothermal circulation 
requires the presence of hot molten rock to locally elevate temperature and drive thermal convection 
and is facilitated by the presence of fractures and faults that allow seawater to penetrate the crust and 
interact with hot rock. Similar to previous work hydrothermal cooling is implemented using an 
effective thermal conductivity that is a function of an effective Nusselt number ``N_{eff}`` and parameters 
that describe the reduction in hydrothermal circulation with depth and temperature due to a reduction
in permeability. The effective Nusselt number ``Nu_{eff}`` is used to increase effective thermal 
conductivity close to areas with molten rock and in regions with relatively high plastic strain rates
and high plastic strain where faults and damaged shear zones are likely to be present. The effective
thermal conductivity model used in this work is given by the following equation:

###### eq:hydrothermal-k
```math
k_{m} = 
\begin{cases}
    k_m^o + k_m^o(Nu_{eff} - 1)
        \exp\left[f_{s}\left(
                        2 - \frac{T_m}{T_{hydro}^{max}} - \frac{ d_{sm,m} }{ d_{mydro}^{max} }
                    \right)
            \right] 
                & \text{if } 0 < d_{sm(m)} < d_{hydro}^{max} \text{ and } y_{ml} > y_{sl} \\
                & \text{ and } y_m > y_{sl} \text{ and } d_{sed,x_m} < 2500m \\
    k_m^o & \text{otherwise}
\end{cases}
```

where ``f_{s}`` is a smoothing factor set equal to ``0.75``, ``T_{hydro}^{max}`` is the maximum temperature
of the hydrothermal system set equal to ``600^{\circ}C``, ``d_{sm,m}`` is the the submud depth of the
marker and ``d_{hydro}^{max}`` is the maximum submud depth of the hydrothermal system set equal to
``4000m``, ``y_{ml}`` is the y-coordinate of the mudline at the x-coordinate of the marker, ``y_{sl}`` is
the y-coordinate of the sea level, ``y_m`` is the y-coordinate of the marker, ``d_{sed,x_m}`` is the
sediment thickness at the x-coordinate of the marker and ``k_m^o`` is the unmodified thermal 
conductivity of the marker. 

The effective Nusselt number ``Nu_{eff}`` is a functions of the distance from the molten zone,
plastic strain rate and plastic strain as described with the following equations:

###### eq:effective-nusselt
```math
\begin{split}
    Nu_{\dot \epsilon} = 
        1 + (Nu_{dist} - 1)
            \left(1 
                - \exp\left(-\frac{\dot \epsilon_{plastic,m}}{\dot \epsilon_{hydro}}\right)
            \right) \text{,}\\
    Nu_{\epsilon} = 
        1 + (Nu_{dist} - 1)
            \left(1 
                - \exp\left(-\frac{\epsilon_{plastic,m}}{\epsilon_{hydro}}\right)
            \right) \text{,}\\
    Nu_{eff} = \max \left( Nu_{\dot \epsilon}, Nu_{\epsilon} \right)
\end{split}
```

where ``\dot \epsilon_{plastic,m}`` is the plastic strain rate of the marker ``m``, 
``\dot \epsilon_{hydro}`` is the hydrothermal reference strain rate set equal to ``10^{-14}s^{-1}``,
``\epsilon_{plastic,m}`` is the plastic strain of the marker ``m``, ``\epsilon_{hydro}`` is the 
hydrothermal plastic strain reference set equal to 0.5, ``Nu_{\dot \epsilon}`` is the Nusselt number
based on the plastic strain rate, ``Nu_{\epsilon}`` is the Nusselt number based on the plastic strain,
and ``Nu_{dist}`` is the distance-based Nusselt number given by:

###### eq:distance-based-nusselt
```math
Nu_{dist} = max \left( 1.0, f_{hydro}Nu_{ref}\right)
```

where ``Nu_{ref}`` is the reference Nusselt number set equal to 2 in this study and ``f_{hydro}`` is
a hydrothermal distance factor given by:

###### eq:hydrothermal-distance-factor
```math
f_{hydro} = \exp \left( - \frac{D_{molten}}{\lambda_{hydro}} \right)
```

where ``\lambda_{hydro}`` is the hydrothermal distance decay length set equal to ``25km``, and ``D_{molten}``
is the distance from the molten zone calculated using the following equation:

###### eq:distance-from-molten-zone
```math
D_{molten} =
\begin{cases}
    0.0 & \text{if } x_{molten}^{min} > x_m > x_{molten}^{max} \\
    x_{molten}^{min} - x_m & \text{if } x_m < x_{molten}^{min} \\
    x_m - x_{molten}^{max} & \text{if } x_m > x_{molten}^{max}
\end{cases}
```

where ``x_m`` is the x-coordinate of the marker and ``x_{molten}^{min}`` and ``x_{molten}^{max}`` are the 
x-limits of the molten zone given by:

###### eq:x-molten-limits
```math
\begin{split}
    x_{molten}^{min} = x_{molten,k} - 0.5W_{molten,k} - D_{molten}^{buffer} \text{,}\\
    x_{molten}^{max} = x_{molten,k} + 0.5W_{molten,k} + D_{molten}^{buffer}
\end{split}
```

where ``W_{molten,k}`` is the width of the gabbroic molten zone in drainage basin ``k``, ``x_{molten,k}`` is 
the x-coordinate of the mid-point of the gabbroic molten zone in drainage basin ``k``, and 
``D_{molten}^{buffer}`` is the buffer distance for the molten zone set equal to ``25km``.