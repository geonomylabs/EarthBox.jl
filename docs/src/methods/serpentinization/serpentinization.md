# Serpentinization

Serpentinization is modeled using the approach of [mezri24](@citet) whereby the serpentinization ratio 
``R_{serp,m}`` is updated as follows:

###### eq:serpentinization-ratio
```math
R_{serp,m} 
    = min\left(
        R_{serp,m}^o + S_{rate}\left(1-R_{serp}^o\right)\Delta t
        \text{,} \quad
        1 - R_{serp}^o
        \right)
```

where ``m`` is the index of a given marker, ``R_{serp,m}^o`` is the previous serpentinization ratio, and 
``S_{rate}`` is the serpentinization rate defined by the following equation:

###### eq:serpentinization-rate
```math
S_{rate} = S_{rate}^{max}f_{serp,T}f_{serp(\epsilon)}
```

where ``S_{rate}^{max}`` is the maximum serpentinization rate, ``f_{serp,T}`` is the temperature factor 
and ``f_{serp(\epsilon)}`` is the plastic strain rate factor defined by the following equations:

###### eq:serpentinization-factors
```math
\begin{split}
    f_{serp,T} =
        \exp\left(
            -\frac{\left(\frac{T_m - T_{serp}^o}{94}\right)^2}{2}
            \right)
        \left(
            1 + \mathrm{erf}
                \left(
                    -4\frac{\left(\frac{T_m - T_{serp}^o}{w_{serp}}\right)}{\sqrt{2}}
                \right)
        \right)
        \quad \text{,} \\
    f_{serp,\dot \epsilon} = 1 - \exp \left(
            -\frac{\dot \epsilon_{plastic,m}}{\dot \epsilon_{serp}}
        \right)
        \quad \text{.}
\end{split}
```

where ``T_{serp}^o`` is the reference serpentinization temperature and ``\dot \epsilon_{serp}``
is the reference serpentinization strain rate. The parameter values used in this work 
for equations [Eq.](@ref eq:serpentinization-rate) and [Eq.](@ref eq:serpentinization-factors) are
provided in [Tab.](@ref tab:serpentinization-example-parameters).

The serpentinization ratio is used to update marker density ``\rho_m`` using the following equation:

###### eq:density-serpentinization
```math
\begin{split}
    f_{\rho,reduction} = 
        1 
        + \left(
            \frac{\left(\rho_{mantle}^{char} - \Delta \rho_{serp}\right)}{\rho_{mantle}^{char}} - 1
        \right)R_{serp,m} \text{,} \\
    \rho_m = f_{\rho,reduction} \rho_m^o
\end{split}
```

where ``\rho_m^o`` is the marker density prior to this update step, ``\rho_{mantle}^{char}`` is the 
characteristic mantle density and ``\Delta \rho_{serp}`` is the density reduction due to serpentinization.

Marker thermal conductivity ``k_m`` is updated for serpentinization using the following equation:

###### eq:thermal-conductivity-serpentinization
```math
k_m = (1 - R_{serp,m})k_m^o  + R_{serp,m}k_{serp}
```

where ``k_m^o`` is the thermal conductivity of the marker prior to the update step and ``k_{serp,m}`` is 
calculated using the following equation:

###### eq:thermal-conductivity-serpentinite
```math
k_{serp} = (
    358.0*(1.0227*k_{serp,20^\circ C} - 1.882)
    *(1.0/T_{K(m)} - 0.00068)
    + 1.84
)
```

where ``k_{serp,20^\circ C}`` is the thermal conductivity of the serpentine at 20``^\circ``C and 
``T_{K(m)}`` is the temperature of marker ``m`` in Kelvin. The parameter values used in this work for 
equations [Eq.](@ref eq:serpentinization-rate), [Eq.](@ref eq:serpentinization-factors) and
[Eq.](@ref eq:thermal-conductivity-serpentinite) are provided in 
[Tab.](@ref tab:serpentinization-example-parameters).

| Parameter | Symbol | Units | Value |
|:----------|:------:|:-----:|------:|
| Maximum serpentinization rate | ``S_{rate}^{max}`` | s``^{-1}`` | 10``^{-11}`` |
| Maximum submud depth of serpentinization | ``d_{serp}^{max}`` | m | 4000 |
| Maximum temperature of serpentinization | ``T_{serp}^{max}`` | ``^{\circ}``C | 340.5 |
| Reference strain rate for serpentinization | ``\dot{\epsilon}_{serp}`` | s``^{-1}`` | 10``^{-13}`` |
| Characteristic mantle density | ``\rho_{mantle}^{char}`` | kg m``^{-3}`` | 3300 |
| Density reduction due to serpentinization | ``\Delta \rho_{serp}`` | kg m``^{-3}`` | 500 |
| Reference thermal conductivity of serpentinite at 20``^\circ``C | ``k_{20^\circ C}`` | W m``^{-1}`` K``^{-1}`` | 2.6 |

#### tab:serpentinization-example-parameters
*Example parameters for serpentinization model.*