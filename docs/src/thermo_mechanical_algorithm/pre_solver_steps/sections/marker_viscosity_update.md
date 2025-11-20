# Marker Creep Viscosity Update Steps

[1] Update non-linear creep viscosity ``\eta_{creep,m}`` for all markers ``m`` 
([Non-linear Viscosity Update Steps](@ref)).

[2] Initialize visco-plastic viscosity ``\eta_{vp,m}`` with creep viscosity ``\eta_{creep,m}``.

[3] Update marker visco-plastic viscosity ``\eta_{vp,m}`` for each marker ``m`` to account for 
extractable melt fraction using equation [eq:viscosity-partial-melt](@ref) with the current 
visco-plastic viscosity as the initial value ``\eta_{vp,m}^o`` and extractable melt fraction 
``M_{extractable,m}``.


## Non-linear Viscosity Update Steps

**For each marker ``m`` do:**

[1] Initialize the deviatoric stress invariant ``\sigma_{II,A,m}'`` using equation 
[eq:second-invariant-stress](@ref) with current marker deviatoric shear stress ``\sigma_{xy,m}'`` and
deviatoric normal stress ``\sigma_{xx,m}'``.

[2] Calculate effective creep viscosity ``\eta_{creep,m}`` using equation 
[eq:composite-viscosity](@ref) with deviatoric stress invariant ``\sigma_{II,A,m}'`` and 
``\dot \epsilon_{eff}`` calculated with current marker temperature ``T_m``, pressure ``P_m``
and deviatoric stress invariant ``\sigma_{II,A,m}'``.
    
[3] Forecast the deviatoric stress invariant ``\sigma_{II,B,m}'`` using equation 
[eq:visco-elastic-stress-forecast](@ref) with current creep viscosity ``\eta_{creep,m}``. 

[4] Update ``\eta_{creep,m}`` using bisection ([Update Effective Viscosity Using Bisection](@ref)).

[5] Copy updated marker creep viscosity ``\eta_{creep,m}`` and copy to 
``\eta_{vp(m)}`` to initialize visco-plastic marker viscosity for time step.


## Update Effective Viscosity Using Bisection

**While ``N_{iterations}`` ``<`` ``N_{max}`` and ``\Delta \sigma > 1`` Do**

[1] Calculate average stress invariant ``\sigma_{II,avg,m}'`` using the following equation:

###### eq:stress-inv-avg-bisection
```math
\sigma_{II,avg,m}' = \frac{1}{2} 
    \left( \sigma_{II,A,m}' + \sigma_{II,B,m}' \right)
```

[2] Calculate effective creep viscosity ``\eta_{creep,m}`` using equation 
[eq:composite-viscosity](@ref) with deviatoric stress invariant ``\sigma_{II,avg,m}'`` and
``\dot \epsilon_{eff}`` calculated with marker temperature ``T_m``, pressure ``P_m`` and 
deviatoric stress invariant ``\sigma_{II,avg,m}'``.

[3] Update the forecasted visco-elastic stress invariant ``\sigma_{II,fc,m}'`` 
using equation [eq:visco-elastic-stress-forecast](@ref) with ``\eta_{creep,m}`` and
current marker deviatoric shear stress ``\sigma_{xy,m}'`` and deviatoric normal stress 
``\sigma_{xx,m}'``.

[4] Update bisection limits:

!!! ukw ""
    ``\text{if } ( \sigma_{II,B,m}' > \sigma_{II,m,A} ) \text{ and } ( \sigma_{II,fc,m}' > \sigma_{II,avg,m}' ) \text{ or } ( \sigma_{II,B,m}' < \sigma_{II,A,m}' ) \text{ and } ( \sigma_{II,fc,m}' < \sigma_{II,avg,m}' )``

    ``\sigma_{II,A,m}' = \sigma_{II,avg,m}'``

    ``\text{else}``

    ``\sigma_{II,B,m}' = \sigma_{II,avg,m}'``

    ``\text{end}``

[5] Update the stress change ``\Delta \sigma`` using the following equation:

###### eq:bisection-stress-change
```math
\Delta \sigma = \sigma_{II,B,m}' - \sigma_{II,A,m}'
```