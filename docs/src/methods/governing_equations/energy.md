# Conservation of Energy

Temperature is obtained by solving the conservation of energy equation as given by:

###### eq:energy
```math
\rho C_p \frac{DT}{Dt}
= \frac{\partial}{\partial x} \left( k \frac{\partial T}{\partial x} \right)
+ \frac{\partial}{\partial y} \left( k \frac{\partial T}{\partial y} \right)
+ H_{rad} + H_{shear} + H_{adi} + H_{melt} + H_{exo}
```

where ``T`` is the temperature, ``C_p`` is the specific heat capacity, ``k`` is the thermal
conductivity, ``H_{rad}`` is the radiogenic heat production, ``H_{shear}`` is the shear heating, 
``H_{adi}`` is the adiabatic heating, ``H_{melt}`` is the latent heat of melting and 
crystallization, ``H_{exo}`` is the heat produced from from exothermic serpentinization reactions 
and ``\frac{DT}{Dt}`` is the substantive time temperature derivative defined as:

###### eq:substantive-temperature-derivative
```math
\frac{DT}{Dt} = \left(\frac{\partial T}{\partial t} + v_x \frac{\partial T}{\partial x} + 
v_y \frac{\partial T}{\partial y}\right)
```

The thermal conductivity term ``k`` in [Eq.](@ref eq:energy) is described with the following 
temperature-dependent model:

###### eq:thermal-conductivity
```math
k = (358(1.0227k_{20^\circ C} - 1.882)(1.0/T_K - 0.00068) + 1.84)
```

where ``k_{20^\circ C}`` is the thermal conductivity of the rock at 20``^\circ``C, ``T_K`` is 
temperature in Kelvin [hantschel09](@citep). The heat capacity term ``C_p`` in 
[Eq.](@ref eq:energy) is described with the following temperature-dependent model:

###### eq:heat-capacity
```math
C_p =
C_{p20}
    \Big(
        0.953 + 2.29 \cdot 10^{-3} T_{^{\circ} C}
        - 2.835 \cdot 10^{-6} \left(T_{^{\circ}C}\right)^2
        + 1.191 \cdot 10^{-9} \left(T_{^{\circ}C}\right)^3
    \Big)
```

where ``C_{p20}`` is the heat capacity of the rock at ``20^\circ C`` and ``T_{^{\circ}C}`` is the 
temperature in ``^{\circ}C`` [hantschel09](@citep).

The shear heating term ``H_{shear}`` is defined using the following equation: 

###### eq:shear-heating
```math
H_{shear} = \frac{\sigma_{xx}'^{2}}{\eta_{vp}} + \frac{\sigma_{xy}'^{2}}{\eta_{vp}}
```

where ``\sigma_{xx}'`` and ``\sigma_{xy}'`` are the deviatoric components of the stress tensor and
``\eta_{vp}`` is the visco-plastic effective viscosity. The adiabatic heating term ``H_{adi}`` is 
defined using the following equation:

###### eq:adiabatic-heating
```math
H_{adi} = c_{adi} \frac{DP}{Dt}
```

where ``\frac{DP}{Dt}`` is the substantive pressure derivative and ``c_{adi}`` is the adiabatic 
coefficient given by:

###### eq:adiabatic-coefficient
```math
c_{adi} = \alpha T
```

where ``\alpha`` is the thermal expansion coefficient. The substantive pressure derivative 
``\frac{DP}{Dt}`` is defined as:

###### eq:substantive-pressure-derivative
```math
\frac{DP}{Dt} = \left(\frac{dP}{dx}v_x + \frac{dP}{dy}vy\right) \text{.}
```

The melt generation term ``H_{melt}`` is divided into two components:

###### eq:melt-generation
```math
H_{melt} = H_{melt, adi} + H_{melt, T}
```

where ``H_{melt, adi}`` is the adiabatic component of the latent heat of melting associated with
changes in pressure within partially molten domains and ``H_{melt, T}`` is the temperature-dependent 
component. The adiabatic component of the latent heat of melting and the temperature-dependent 
components are given by:

###### eq:melt-generation-adiabatic
```math
H_{melt, adi} = -\rho L\frac{\partial M}{\partial P}\frac{DP}{Dt}
```

###### eq:melt-generation-temperature
```math
H_{melt, T} = -\rho L \frac{\partial M}{\partial T}\frac{DT}{Dt}
```

where ``L`` is the latent heat of melting, ``M`` is the melt fraction, and ``\frac{DT}{Dt}`` is the
substantive temperature derivative. Equation [Eq.](@ref eq:melt-generation-temperature) is included 
in an effective heat capacity term in the energy equation [Eq.](@ref eq:energy) to account for the 
temperature-dependent latent heat of melting and crystallization. The effective heat capacity term 
is given by:

###### eq:effective-heat-capacity
```math
C_p^{eff} = C_p + L \frac{\partial M}{\partial T}
```

The partial derivative ``\frac{\partial M}{\partial T}`` in 
[Eq.](@ref eq:effective-heat-capacity) is calculated using a finite difference approximation as 
follows:

###### eq:partial-melt-fraction-temperature
```math
\frac{\partial M}{\partial T} = \frac{M_{T + \Delta T,P} - M_{T - \Delta T,P}}{2\Delta T}
```

where ``\Delta T = 1 K`` and ``M_{T + \Delta T,P}`` and ``M_{T - \Delta T,P}`` are the melt 
fractions at the current pressure ``P`` and temperatures ``T + \Delta T`` and ``T - \Delta T``, 
respectively.

Similarly, equation [Eq.](@ref eq:melt-generation-adiabatic) is included in an effective adiabatic
coefficient term in equation [Eq.](@ref eq:adiabatic-heating) to account for the adiabatic 
component of the latent heat of melting and crystallization. The effective adiabatic coefficient 
term is given by:

###### eq:effective-adiabatic-coefficient
```math
c_{adi}^{eff} = c_{adi} - \rho L \frac{\partial M}{\partial P}
```

The partial derivative ``\frac{\partial M}{\partial P}`` in 
[Eq.](@ref eq:effective-adiabatic-coefficient) is calculated using a finite difference 
approximation as follows:

###### eq:partial-melt-fraction-pressure
```math
\frac{\partial M}{\partial P} = \frac{M_{T,P + \Delta P} - M_{T,P - \Delta P}}{2\Delta P}
```

where ``\Delta P = 1000 Pa`` and ``M_{T,P + \Delta P}`` and ``M_{T,P - \Delta P}`` are the melt 
fractions at the current temperature ``T`` and pressures ``P + \Delta P`` and ``P - \Delta P``, 
respectively. The exothermic heat production term ``H_{exo}`` associated with serpentinization is 
given by:

###### eq:exothermic-heat-production
```math
H_{exo} = \frac{f_{serp,inc}E}{\Delta t M_{serp}}
```

where ``f_{serp,inc}`` is the incremental of serpentinization ratio, ``M_{serp}`` is the molar 
volume of serpentine (``m^3/mol``), ``E`` is the enthalpy change (``J/mol``) and ``\Delta t`` is 
the time step.

The substantive time pressure derivative that appears in [Eq.](@ref eq:adiabatic-heating) and 
[Eq.](@ref eq:melt-generation-adiabatic) is approximated using the following equation that 
neglects deviations of dynamic pressure gradients from lithostatic pressure gradients 
[gerya2019](@citep):

###### eq:pressure-derivative
```math
\left(\frac{dP}{dx}v_x + \frac{dP}{dy}v_y\right) = \rho g_x v_x + \rho g_y v_y
```