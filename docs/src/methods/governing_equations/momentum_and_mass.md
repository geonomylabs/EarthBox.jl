# Conservation of Momentum and Mass

Velocity and pressure are obtained by solving the conservation momentum and mass equations for a
slow, highly viscous incompressible visco-elasto-plastic fluid:

###### eq:x-stokes
```math
\frac{\partial \sigma_{xx}'}{\partial x}
+ \frac{\partial \sigma_{xy}'}{\partial y}
- \frac{\partial P}{\partial x}
= -\rho g_x
```

###### eq:y-stokes
```math
-\frac{\partial \sigma_{xx}'}{\partial y} 
+ \frac{\partial \sigma_{xy}'}{\partial x}
- \frac{\partial P}{\partial y}
= -\rho g_y
```

###### eq:visco-elastic-term
```math
f_{ve} = \frac{\eta_{vp}}{\mu \Delta t + \eta_{vp}}
```

###### eq:xx-stress
```math
\sigma_{xx}'
= 2 \eta_{vp} \dot{\epsilon}_{xx}'(1 - f_{ve}) + \sigma_{xx}'^{co} f_{ve}
```

###### eq:xy-stress
```math
\sigma_{xy}'
= 2 \eta_{vp} \dot{\epsilon}_{xy} (1 - f_{ve}) + \sigma_{xy}'^{co} f_{ve}
```

###### eq:strain-rate-components
```math
\dot{\epsilon}_{xx} = \frac{\partial v_x}{\partial x}, \quad 
\dot{\epsilon}_{xy} = \frac{1}{2} \left( \frac{\partial v_x}{\partial y} + 
\frac{\partial v_y}{\partial x} \right)
```

###### eq:continuity
```math
\frac{\partial v_x}{\partial x} + \frac{\partial v_y}{\partial y} = 0
```

where ``\sigma_{xx}'`` is the deviatoric normal stresses, ``\sigma_{xx}'^{co}`` is the corrected
deviatoric normal stress from the previous time step, ``\sigma_{xy}'`` is the deviatoric 
shear stress, ``\sigma_{xy}'^{co}`` is the corrected deviatoric shear stress from the previous time 
step, ``P`` is pressure, ``\rho`` is density, ``g_x`` and ``g_y`` are the gravitational 
accelerations, ``\eta_{vp}`` is the visco-plastic viscosity, ``\mu`` is the shear modulus, 
``\Delta t`` is the computational time step, ``\dot{\epsilon}_{xx}'`` is normal deviatoric strain 
rate, ``\dot{\epsilon}_{xy}'`` is the deviatoric shear strain rate and ``v_x`` and ``v_y`` are the 
velocity components. Equations [Eq.](@ref eq:x-stokes) and [Eq.](@ref eq:y-stokes) are the x- and 
y-components of the Stokes equations, respectively. We note that the equation [Eq.](@ref eq:y-stokes) 
utilizes the relationships ``\sigma_{yy}' = -\sigma_{xx}`` and ``\sigma_{yx}' = \sigma_{xy}`` to 
reduce computational storage requirements. The Boussinesq approximation is used to account for the 
effects of temperature and pressure on density in the gravitation terms in [Eq.](@ref eq:x-stokes) and 
[Eq.](@ref eq:y-stokes) while maintaining the assumptions of incompressibility from 
[Eq.](@ref eq:continuity):

###### eq:density
```math
\rho = \rho_{ref} \left(1 - \alpha (T - T_{ref})\right)\left(1 + \beta(P-P_{ref})\right)
```

where ``\rho_{ref}`` is the reference density, ``\alpha`` is the thermal expansion coefficient,
``T_{ref}`` is the reference temperature (298.15 K), ``\beta`` is the compressibility coefficient,
``P_{ref}`` is the reference pressure (``10^5 Pa``). Gravitational acceleration in the x-direction is set 
equal to zero, ``g_x = 0``, and in the y-direction gravitational acceleration is set to ``g_y = 9.8`` 
m/s``^2``. The visco-plastic viscosity term ``\eta_{vp}`` is given by:

###### eq:visco-plastic-viscosity
```math
\begin{split}
\eta_{vp} = \eta_{creep} \quad \text{for} \quad \sigma_{II, elastic} \leq \sigma_{yield} 
\text{, and} \\
\eta_{vp} = \mu \Delta t \frac{\sigma_{yield}}{\sigma_{II, elastic} - \sigma_{yield}} \quad 
\text{for} \quad \sigma_{II, elastic} > \sigma_{yield}
\end{split}
```

where ``\eta_{creep}`` is the viscosity associated with creep mechanisms, ``\sigma_{yield}`` 
is the yield stress and ``\sigma_{II, elastic}`` is the second invariant for purely elastic stress
buildup as given by

###### eq:second-invariant-elastic-stress
```math
\sigma_{II, elastic} = \frac{\mu \Delta t + \eta_{creep}}{\eta_{creep}}\sigma_{II}
```

where ``\sigma_{II}`` is the second invariant of the stress tensor given by

###### eq:second-invariant-stress
```math
\sigma_{II} = \sqrt{\left(\sigma_{xy}^2 + \sigma_{xx}^2\right)} \text{.}
```

Visco-plastic viscosity ``\eta_{vp}`` is is limited by a minimum value and a maximum value to
avoid numerical issues. A typical minimum viscosity used in EarthBox simulations is ``10^18 Pa \cdot s``,
and a typical maximum viscosity is ``10^26 Pa \cdot s``.

The yield stress ``\sigma_{yield}`` is given by the following equations:

###### eq:yield-stress
```math
\sigma_{yield} =
\begin{cases}
    \sigma_c \cos{\theta} + \sin{\theta}\left(P - P_f\right) & \quad \text{for} P \geq P_f\\
    \sigma_c \cos(\theta) & \quad \text{for} P < P_f
\end{cases}
```

where ``\sigma_c`` is cohesion, ``\theta`` is the friction angle, and ``P_f`` is 
fluid pressure. Yield stress is limited by a minimum value of 1 MPa.

The corrected deviatoric stresses from previous time steps, ``\sigma_{xx}'^{co}`` and 
``\sigma_{xy}'^{co}``, which appear in equations [Eq.](@ref eq:xx-stress) and [Eq.](@ref eq:xy-stress), 
take into account rotation using the following equations:

###### xx-stress-rot
```math
\sigma_{xx}'^{co}
= \sigma_{xx}'\left(\cos(\omega \Delta t)^{2} - \sin(\omega \Delta t)^{2}\right)
- \sigma_{xy}' \sin(2\omega \Delta t)
```

###### xy-stress-rot
```math
\sigma_{xy}'^{co} = \sigma_{xx}' \sin(2\omega \Delta t) + \sigma_{xy}' \cos(2\omega \Delta t)
```

```math
\omega = \frac{1}{2}\left(\frac{\partial v_y}{\partial x} - \frac{\partial v_x}{\partial y}\right)
```

where ``\omega`` is the rotation rate or spin of the material.

A free surface is implemented in the model using a sticky-air-water layer requiring stabilization to 
avoid the so called drunken sailor effect. To achieve this we use the approach described by 
[gerya2019](@citet) whereby advection-related density changes are incorporated into [Eq.](@ref eq:y-stokes) 
to better approximate density at the end of the time step and avoid the drunken sailor effect. The 
modified y-Stokes equation is given by:

###### eq:y-stokes-stabilization
```math
-\frac{\partial \sigma_{xx}'}{\partial y}
+ \frac{\partial \sigma_{xy}'}{\partial x}
- \frac{\partial P}{\partial y} 
- g_y \Delta t \left( v_x \frac{\partial \rho}{\partial x} + 
v_y \frac{\partial \rho}{\partial y} \right)
= -\rho g_y
```
