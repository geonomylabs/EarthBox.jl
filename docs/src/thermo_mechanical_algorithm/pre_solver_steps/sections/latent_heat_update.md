# Latent Heat Update Steps

[1] Update the adiabatic coefficient ``c_{adi,m}`` to include latent heat using equation
[eq:effective-adiabatic-coefficient](@ref) with current adiabatic coefficient ``c_{adi,m}`` to define
the adiabatic coefficient prior to this update step, latent heat of fusion ``L_m`` and equation
[eq:partial-melt-fraction-pressure](@ref) to calculate the partial derivative of melt fraction
with respect to pressure with marker temperature ``T_m`` and pressure ``P_m``.

[2] Update marker heat capacity ``C_{p,m}`` to include latent heat using equation
[eq:effective-heat-capacity](@ref) with current heat capacity ``C_{p,m}`` to define the heat capacity
prior to this update step, latent heat of fusion ``L_m`` and equation
[eq:partial-melt-fraction-temperature](@ref) to calculate the partial derivative of melt fraction
with respect to temperature with marker temperature ``T_m`` and pressure ``P_m``.