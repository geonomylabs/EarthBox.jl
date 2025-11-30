# Viscous Creep

A composite rheology is used to model solid-state creep including diffusion, dislocation, and Peierls 
creep. The composite effective strain rate invariant ``\dot \epsilon_{eff}`` is given by:

###### eq:composite-strain-rate
```math
\dot \epsilon_{eff} = 
    \max\left(\dot \epsilon_{dis}, \dot \epsilon_{pei}\right) + \dot \epsilon_{dif} \text{. }
```

where ``\dot \epsilon_{dis}``, ``\dot \epsilon_{pei}``, and ``\dot \epsilon_{dif}`` are the strain rate
invariants for dislocation, Peierls, and diffusion creep, respectively. The effective creep 
viscosity ``\eta_{creep,m}`` is then given by:

###### eq:composite-viscosity
```math
\eta_{creep,m} = \frac{\sigma_{II(m)}'}{2\dot \epsilon_{eff}}
```

where ``\sigma_{II(m)}'`` is the second invariant of the stress tensor in ``MPa``.

The strain rate invariant for diffusion creep ``\dot \epsilon_{dif}`` is calculated using the 
following equation:

###### eq:diffusion_creep
```math
\dot \epsilon_{dif} = 
    A_{dif}\exp\left[
        -\frac{
            \left( E_{dif(m)}10^3 + V_{dif(m)}P_m10^{-6} \right)
        }{R_{gas}T_m}
        \right](\sigma_{II(m)}'10^{-6})
```

where ``A_{dif(m)}`` is the diffusion creep pre-exponential factor in ``1/s/MPa``, ``E_{dif(m)}`` is the 
diffusion creep activation energy in ``kJ/mol``, ``V_{dif(m)}`` is the diffusion creep activation 
volume in ``J/MPa/mol``, ``R_{gas}`` is the ideal gas constant in ``J/mol/K``, ``T_m`` is the marker
temperature in ``K``, and ``P_m`` is the marker pressure in ``Pa``. The strain rate invariant for 
dislocation creep ``\dot \epsilon_{dis}`` is given by:

###### eq:dislocation_creep
```math
\dot \epsilon_{dis} = 
    A_{dis}\exp\left[
        -\frac{
            \left( E_{dis(m)}10^3 + V_{dis(m)}P_m10^{-6} \right)
        }{R_{gas}T_m}
        \right](\sigma_{II(m)}'10^{-6})^{n_m}
```

where ``A_{dis}`` is the dislocation creep pre-exponential factor in ``1/s/(MPa)^n``, ``E_{dis}`` is the 
dislocation creep activation energy in ``kJ/mol``, ``V_{dif}`` is the dislocation creep activation 
volume in ``J/MPa/mol`` and ``n_m`` is the stress exponent. finally, the strain rate invariant for 
Peierls creep ``\dot \epsilon_{pei}`` is calculated using the following equation:

###### eq:peierls_creep
```math
\dot \epsilon_{pei} = 
    A_{pei} (\sigma_{II(m)}'10^{-6})^{2} \\
        \exp\left[
        -\frac{\left( E_{pei(m)}10^3 + V_{pei(m)}P_m10^{-6} \right)}{R_{gas}T_m}
        \left(1 - \left( \frac{\sigma_{II(m)}'10^{-6}}{\sigma_{pei(m)}} \right)^{n1_m}\right)^{n2_m}
        \right]
```

where ``A_{pei(m)}`` is the Peierls creep pre-exponential factor in ``1/s/(MPa)^2``, ``E_{pei(m)}`` is the
Peierls creep activation energy in ``kJ/mol``, ``V_{pei(m)}`` is the Peierls creep activation volume in
``J/MPa/mol``, ``\sigma_{pei(m)}`` is the Peierls stress in ``MPa`` and ``n1_m`` and ``n2_m`` are the 
Peierls stress exponents.

The visco-elastic nature of the rheology implemented in this work and the power law and exponential
nature of the creep mechanisms introduces a non-linear time-dependent relationship between stress and 
effective viscosity. Similar to [gerya2010](@cite) a bisection algorithm 
([Viscosity Using Bisection](@ref "Update Effective Viscosity Using Bisection")) is used to solve for the
effective viscosity consistent with a visco-elastic stress forecast for a given marker ``m`` 
calculated using the following equations:

###### eq:visco-elastic-stress-forecast
```math
    f_{ve,m} = \frac{ \eta_{creep,m} }{ \mu_m \Delta t + \eta_{creep,m} }
        \text{,} \\
    \sigma_{xx,fc,m}' = 
        2 \eta_{creep,m} \dot{\epsilon}_{xx,m}'R_{sr,m}(1-f_{ve,m})
        + \sigma_{xx,m}' f_{ve,m} 
                \text{,} \\
    \sigma_{xy,fc,m}' =
        2 \eta_{creep,m} \dot{\epsilon}_{xy,m}'R_{sr,m}(1-f_{ve,m})
        + \sigma_{xy,m}' f_{ve,m}
                \text{,} \\
    \sigma_{II,fc,m}' = max 
        \left(
        \sqrt{
        \left( \sigma_{xx,fc,m}' \right)^2 
        + \left( \sigma_{xy,fc,m}' \right)^2
        }, \sigma_{min}
        \right)
```

where ``\sigma_{min}`` is a user-defined minimum stress value, ``f_{ve}`` is the visco-elastic factor, 
``\mu_m`` is the shear modulus, ``\Delta t`` is the time step and ``\dot{\epsilon}_{xx,m}'`` and 
``\dot{\epsilon}_{xy,m}'`` are the deviatoric strain rate components interpolated from strain rate 
defined on the staggered grid as described by equation [eq:deviatoric-strain-rate](@ref),
and ``R_{sr,m}`` is strain rate ratio given by the following equation:

###### eq:strain_rate_ratio
```math
R_{sr,m} = \frac{\dot{\epsilon}_{II,stress,m}'}{\dot{\epsilon}_{II, velocity,m}'}
```

where ``\dot{\epsilon}_{II,velocity,m}'`` is the strain rate invariant for marker ``m`` calculated using
the following equation:

###### eq:marker_strain_rate_invariant
```math
\dot{\epsilon}_{II,velocity,m}' = \sqrt{
    \left( \dot{\epsilon}_{xx,m}' \right)^2 + \left( \dot{\epsilon}_{xy,m}' \right)^2
} \text{.}
```

and ``\dot{\epsilon}_{II,stress,m}'`` is the strain rate invariant for marker ``m`` calculated using nodal 
deviatoric stress and deviatoric stress changes interpolated from the staggered grid as follows:

###### eq:stress_based_strain_rate
```math
\dot{\epsilon}_{II,stress,m}' = \sqrt{
    \left( \frac{\sigma_{xx,m_2}'}{2\eta_{vp,m}} + \frac{\Delta \sigma_{xx,m}'}{2\Delta t \mu_m} \right)^2
    + \left( \frac{\sigma_{xy,m_2}'}{2\eta_{vp,m}} + \frac{\Delta \sigma_{xy,m}'}{2\Delta t \mu_m} \right)^2
}
```

where ``\sigma_{xx,m_2}'`` and ``\sigma_{xy,m_2}'`` are the forecasted deviatoric stress components at
markers interpolated from forecasted deviatoric grid stress ``\sigma_{xx(i,j)_{p2}}'`` and
``\sigma_{xy(i,j)_{b2}}'`` from equation [eq:deviatoric-stress-forecast](@ref) and ``\Delta \sigma_{xx,m}'`` and 
``\Delta \sigma_{xy,m}'`` are the grid stress changes at markers interpolated from 
``\Delta \sigma_{xx(i,j)_p}'`` and ``\Delta \sigma_{xy(i,j)_b}'`` from equation 
[eq:visco-elasto-plastic-stress-changes](@ref). The use of ``R_{sr,m}`` in equation 
eq:visco-elastic-stress-forecast makes the visco-elastic stress forecast dependent on nodal 
stress and reduces numerical diffusion in shear zones [gerya2010](@citep).

With each time step, the effective visco-plastic viscosity ``\eta_{vp,m}`` is 
initialized with the effective creep viscosity ``\eta_{creep,m}`` before the effects of plastic 
failure and partial melt are taken into account using equation [eq:visco-plastic-viscosity](@ref) 
and the following equation, respectively:

###### eq:viscosity-partial-melt
```math
    \eta_{vp,m} = \eta_{vp,m}^o /\exp\left(f_{\alpha}M_{extractable,m}\right) 
        \quad for \quad M_{extractable,m} < 0.3 \quad \text{and,} \\
    \eta_{vp,m} = \eta_{melt} \quad for \quad M_m \geq 0.3
```

where ``\eta_{vp,m}^o`` is the previous visco-plastic viscosity initialized as ``\eta_{creep,m}``, 
``f_{\alpha}`` is a user-defined parameter that controls the sensitivity of the visco-plastic viscosity to 
melt fraction and ``\eta_{melt}`` is the viscosity of melt at ``M_{mextractable,m} = 1``. The parameter 
``f_{\alpha}`` is set to 35 by default.