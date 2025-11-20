# Pre-solver Steps

[1] Update boundary conditions for velocity stepping and transient temperature along
the lower boundary of the model.

[2] Randomize initial marker friction angles ``\theta_m^o`` using equation 
[eq:friction-angle-randomization](@ref).

[3] Update marker composition ``C_m`` to account for compositional effects associated 
with fractionation of normal gabbroic melt by transforming gabbroic magma markers to 
fractionated gabbroic melt if the distance between the marker and the Moho is 2km or 
less.

[4] Update marker melt fraction ``M_m`` based on marker temperature ``T_m`` and marker 
pressure ``P_m`` using equation [eq:melt-fraction](@ref).

[5] Update marker composition ``C_m`` by transforming purely molten markers that have 
cooled below the liquidus temperature to a solid state using equations 
[eq:liquidus-mantle](@ref), [eq:solidus-liquidus-gabbro](@ref) and 
[eq:solidus-liquidus-fractionated-gabbro](@ref).

[6] Update marker composition ``C_m`` by transforming solid markers with temperature 
between the solidus and liquidus temperatures to a partially molten state using 
equations [eq:solidus-mantle](@ref), [eq:liquidus-mantle](@ref), 
[eq:solidus-liquidus-gabbro](@ref) and [eq:solidus-liquidus-fractionated-gabbro](@ref).

[7] Update marker composition ``C_m`` by transforming partially molten markers that have 
cooled below the solidus to a solid state by updating ``C_m`` using equations 
[eq:solidus-mantle](@ref), [eq:solidus-liquidus-gabbro](@ref) and 
[eq:solidus-liquidus-fractionated-gabbro](@ref).

[8] Update rock properties for density ``\rho_m``, heat capacity ``C_{p,m}``, thermal 
conductivity ``k_m``, adiabatic coefficient ``c_{adi,m}`` for each marker ``m`` taking into 
account temperature, pressure, sediment porosity, serpentinization, partial melt and 
melt depletion ([Rock Property Update Steps](@ref)).

[9] Include the effects of latent heat in the energy equation by adding latent heat 
terms in the adiabatic coefficient ``c_{adi,m}`` and heat capacity ``C_{p,m}`` for each 
marker ``m`` ([Latent Heat Update Steps](@ref)).

[10] Include the effects of hydrothermal cooling on regional thermal structure by 
modifying the marker thermal conductivity ``k_m`` 
([Hydrothermal Update Steps](@ref)).

[11] Initialize bilinear interpolation parameters for each marker ``m`` including indices 
and distances for the upper-left basic node and weights of the four basic nodes of 
the cell containing marker ``m`` and summed marker weights on basic nodes, pressure nodes 
and ``v_y`` grid nodes used in the denominator of equation 
[eq:bilinear_interp2grid](@ref).

[12] Interpolate thermal parameters from markers to basic grids:

```math
\begin{split}
\text{Temperature (inclusive): } & T_m \rightarrow T_{(i,j)_b} \\
\text{Density x Heat Capacity (inclusive): } & \rho_m C_{p,m} \rightarrow 
\rho_{(i,j)_b} C_{p(i,j)_b} \\
\text{Thermal conductivity (inclusive): } & k_m \rightarrow k_{(i,j)_b} \\
\text{Radiogenic Heat Production (inclusive): } & H_m^{rad} \rightarrow 
H_{(i,j)_b}^{rad} \\
\text{Adiabatic Coefficient (inclusive): } & c_{adi,m} \rightarrow c_{adi(i,j)_b}
\end{split}
```

[13] Update marker failure properties including cohesion ``\sigma_{c,m}`` and friction 
angle ``\theta_m`` including the effects of randomized initial friction angle ``\theta_m^o``, 
strain weakening and melt damage ([Marker Failure Properties Update Steps](@ref)).

[14] Update marker viscosity including non-linear effective creep viscosity 
``\eta_{creep,m}`` and initialized visco-plastic viscosity ``\eta_{vp,m}`` taking into 
account partial melt ([Marker Creep Viscosity Update Steps](@ref)).

[15] Update melt damage factor ``f_{damage(m)}`` ([Melt Damage Update Steps](@ref)).

[16] Interpolate Stokes-continuity parameters from markers to staggered grids:

```math
\text{Density (inclusive): } \rho_m \rightarrow \rho_{(i,j)_b} \\
\text{Cohesion (exclusive): } C_m \rightarrow C_{(i,j)_b} \\
\text{Friction Angle (exclusive): } \theta_m \rightarrow \theta_{(i,j)_b} \\
\text{Creep Viscosity (exclusive): } \eta_{creep,m} \rightarrow \eta_{creep{(i,j)_b}} \\
\text{Visco-Plastic Shear Viscosity (exclusive): } \eta_{vp,m} \rightarrow \eta_{vp(i,j)_b} \\
\text{Visco-plastic Normal Viscosity (exclusive): } \eta_{vp,m} \rightarrow \eta_{vp(i,j)_p} \\
\text{Shear Modulus on basic nodes (exclusive): } \mu_m \rightarrow \mu_{(i,j)_b} \\
\text{Shear modulus on pressure nodes (inclusive): } \mu_m \rightarrow \mu_{(i,j)_p} \\
\text{Deviatoric shear stress (exclusive): } \sigma_{xy,m}' \rightarrow \sigma_{xy(i,j)_{b1}}' \\
\text{Deviatoric normal stress (inclusive): } \sigma_{xx,m}' \rightarrow \sigma_{xx(i,j)_{p1}}' \\
```