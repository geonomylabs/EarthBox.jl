# Rock Property Update Steps

[1] Update density ``\rho_m`` for each marker ``m`` using equation [eq:density](@ref)
with thermal expansion coefficient ``\alpha_m``, compressibility ``\beta_m``, temperature ``T_m``, 
pressure ``P_m``, and reference density ``\rho_{ref,m}``.

[2] Update heat capacity ``C_{p,m}`` for each marker ``m`` using equation [eq:heat-capacity](@ref) 
with heat capacity at ``20^\circ C`` ``C_{p20,m}`` and marker temperature in degrees Celsius 
``T_{^{\circ}C,m}``.

[3] Update marker thermal conductivity ``k_m`` for each marker ``m`` using equation 
[eq:thermal-conductivity](@ref) with reference thermal conductivity ``k_{20^\circ C,m}`` and
temperature ``T_m``.

[4] Update marker adiabatic coefficient ``c_{adi,m}`` using equation 
[eq:adiabatic-coefficient](@ref) with thermal expansion coefficient ``\alpha_{m}`` and marker 
temperature ``T_m``.

[5] Update sediment rock properties ``\rho_m``, ``k_m`` and ``C_{p,m}`` for each marker ``m`` to account 
for sediment porosity ``\phi_m`` if marker composition ``C_m`` is sediment using equations 
[eq:sediment-rock-props](@ref) with matrix properties for density, thermal conductivity and heat 
capacity equal to the current marker properties ``\rho_m``, ``k_m`` and ``C_{p,m}``, respectively.

[6] Update ``\rho_m`` and ``k_m`` for each marker ``m`` to account for serpentinization using equations
[eq:density-serpentinization](@ref) and [eq:thermal-conductivity-serpentinization](@ref) with
serpentinization ratio ``R_{serp,m}`` and current marker properties ``\rho_m`` and ``k_m`` to define
properties prior to update step.

[7] Update marker density ``\rho_m`` for each marker ``m`` to account for depletion of rock matrix 
due to partial melting using equation [eq:density-melt-depletion](@ref) with extracted melt 
fraction ``M_{extracted,m}`` and current marker density ``\rho_m`` used to define density prior to 
this update step.

[8] Update marker density ``\rho_m`` for each marker ``m`` to account for partial melt in pore 
space using equation [eq:density-melt-pore-space](@ref) with pressure ``P_m``, extractable melt
fraction ``M_{extractable,m}`` and the current marker density ``\rho_m`` to define density prior to
this update step.
