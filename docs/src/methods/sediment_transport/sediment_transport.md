# Sediment Transport

Sediment transport is modeled on a Eulerian topography grid ``(x_{topo,i}, y_{topo,i})`` where ``i`` ranges
from ``0`` to ``N_t-1``, ``N_{topo}`` is the number of topography nodes and ``y_{t,i}`` is the y-coordinate of 
the topography at the ``i``-th node. The spacing of the topography grid in the x-direction is denoted by 
``\Delta x_t``, which remains constant throughout the model domain. Prior to applying the sediment 
transport model, topography grid nodes are advected in the current velocity field using a 4th order 
Runge-Kutta scheme and then interpolated to the original topography grid to obtain new y-coordinates 
of topography that track the interface between rock and the sticky-air/water layer.

Sediment transport is modeled using a diffusion equation that assumes density does not differ 
between sediment and bedrock and that the cohesion of bedrock is negligible as to not limit sediment
supply at the surface:

###### eq:sediment-transport
```math
\frac{\partial y_{topo}}{\partial t} = \frac{\partial}{\partial x}
    \left( \kappa_s \frac{\partial y_{topo}}{\partial x} \right) + R_{pelagic}
```

where ``y_{topo}`` is the y-coordinate of topography in the 2D model domain, ``\kappa_s`` is the 
sediment transport diffusivity, and ``R_{pelagic}`` is a sediment source term that accounts for pelagic 
sedimentation and sediment transport in and out of the plain of the model domain. We use the approach 
described by [martinez19](@citet) to define ``\kappa_s`` so that it includes the effects of subaerial 
slope diffusion, subaerial fluvial concentrative diffusional transport associated with river systems 
within a given drainage basin and water-depth dependent submarine slope diffusion that accounts for 
wave and tidal effects at shallow water depths:

###### eq:sediment-transport-diffusivity
```math
\kappa_{s,i} =
\begin{cases}
    \kappa_{sm} \exp\left( -\frac{W}{\lambda_{sm}} \right)
        & \text{if } W_i > 0 \\
    \kappa_{sa} + \kappa_{fluvial,i}
        & \text{otherwise}
\end{cases}
```

where ``W_i`` is water depth in ``m``, ``\kappa_{sm}`` is the submarine slope diffusivity in ``m^2/s``, 
``\lambda_{sm}`` is the submarine diffusivity water-depth decay term in ``m``, ``\kappa_{sa}`` 
is the subaerial slope diffusivity in ``m^2/s`` and ``\kappa_{fluvial}`` is the diffusivity associated with
fluvial transport processes in ``m^2/s``. The fluvial diffusivity ``\kappa_{fluvial,i}`` is described by 
the following equation:

```math
\kappa_{fluvial,i} = R_{precip}C_{subaerial}D_{stream,i}
```

where ``R_{precip}`` is the precipitation rate in ``m/s``, ``C_{subaerial}`` is the subaerial transport coefficient,
and ``D_{stream,i}`` is the downstream distance in ``m`` between drainage divides to the left and right of node ``i``. 
The downstream distance ``D_{stream,i}`` essentially scales fluvial diffusivity with the size of the drainage 
basin. A 1D conservative finite element method is used to solve the sediment transport equation 
[Eq.](@ref eq:sediment-transport) with a time step of ``\Delta t_{sed}``.

The sediment transport equation [Eq.](@ref eq:sediment-transport) does not take into account the effects of
sediment compaction and assumes that sediment density is equal to bedrock density. This assumption
is equivalent to assuming that the sediment is deposited in a state of maximum compaction and that
the substrate upon which new sediment is deposited does not compact. We apply a compaction correction 
to the topography solution at each time step of the sediment transport solver to account for deposition
of sediment in a less compaction state and the compacted of pre-existing sediment:

###### eq:sediment-compaction-correction
```math
y_{topo,i} = 
    y_{topo,trans,i} 
    - \left[
        \left( H_{sed,i}^f + \Delta H_{sed,i}^f\right)
        - \left( H_{sed,i}^o + \Delta H_{sed,i}^o \right) 
    \right]
```

where ``y_{topo,trans,i}`` is the topography obtained by solving the sediment transport equation
[Eq.](@ref eq:sediment-transport) at node ``i``, ``H_{sed,i}^o`` is the initial sediment thickness prior to 
the deposition of new sediment, ``\Delta H_{sed,i}^o`` is the newly deposited sediment thickness in 
the zero-porosity maximum compaction state, ``\Delta H_{sed,i}^f`` is the newly deposited sediment 
thickness in a de-compacted state and ``H_{sed,i}^f`` is the initial sediment thickness compacted to a 
depth of ``\Delta H_{sed,i}^f``. The parameter ``\Delta H_{sed,i}^o`` is related to the new and old 
topography solutions obtained from solving [Eq.](@ref eq:sediment-transport) as follows:

###### eq:depo-thickness-max-compaction
```math
\Delta H_{sed,i}^o = \max\left[ (y_{topo,i}^o - y_{topo,trans,i}) \text{, } 0 \right] \text{.}
```

The newly deposited decompacted sediment thickness ``\Delta H_{sed,i}^f`` is calculated by de-compacting 
``\Delta H_{sed,i}^o`` from an assumed depth below mudline of 12 km to 0 km using the 
conservation of mass and assuming that the relationship between porosity and depth in the sediment 
column can be described by the following equation:

###### eq:sediment-porosity
```math
\phi = \phi_o \exp \left(-\frac{y_{sm,max}}{\lambda_{comp}} \right)
```

where ``\phi`` is the sediment porosity, ``\phi_o`` is the initial sediment porosity at the mudline, 
``y_{sm,max}`` is the maximum depth below the seafloor encountered by a particle during burial 
and ``\lambda_{comp}`` is the compaction decay length [sclater80](@citep). A similar approach is used to
obtain ``H_{sed,i}^f`` whereby equation [Eq.](@ref eq:sediment-porosity) is used to compact ``H_{sed,i}^o`` to 
a depth below mudline of ``\Delta H_{sed,i}^f``. With each time step of the sediment transport solver a correction is
also applied to the location of sediment and sticky air/water markers whereby pre-existing sediment 
and sticky air/water markers are advected vertically using the compaction displacement field.

The marker composition field is updated to account for erosion and sedimentation processes after
the sediment transport solver and compaction correction have been applied to update the y-coordinate 
of topography ``y_{topo,i}``. The updated topography is first interpolated at the x-coordinate ``x_m`` of 
marker ``m`` to determine the elevation of the topography at the marker ``y_{topo,m}``. If the 
y-coordinate of the marker ``y_m`` is less than ``y_{topo,m}`` and the marker is lithological, the 
marker is transformed to either a sticky-air or sticky-water marker depending on whether the marker 
is above or below sea level to account for erosion. If the marker is sticky-air or sticky-water and 
``y_m`` is greater than ``y_{topo,m}``, the marker is transformed to a sedimentary marker.

For a given marker ``m``, sediment density ``\rho_m``, thermal conductivity ``k_m`` and heat capacity ``C_{p,m}`` 
are updated to account for water filled porosity as described in equation [Eq.](@ref eq:sediment-porosity) 
using the following equations:

###### eq:sediment-rock-props
```math
\begin{split}
    \rho_m = (1 - \phi_m) \rho_{matrix,m} + \phi_m \rho_{water} \text{, } \\
    k_m = (1 - \phi_m) k_{matrix,m} + \phi_m k_{water} \\
    C_{p,m} = (1 - \phi_m) C_{p,matrix,m} + \phi_m C_{p,water}
\end{split}
```

where ``\rho_{matrix,p}`` is the density of the rock matrix, ``k_{matrix,m}`` is the thermal conductivity of 
the rock matrix, ``C_{p,matrix,m}`` is the heat capacity of the rock matrix, ``\rho_{water}`` is the density 
of water, ``k_{water}`` is the thermal conductivity of water and ``C_{p,water}`` is the heat capacity of water.

| Parameter | Symbol | Units | Value |
|:----------|:------:|:-----:|------:|
| Subaerial Transport coefficient | ``C_{subaerial}`` | - | ``10^{-4}`` |
| Subaerial slope diffusivity | ``\kappa_{submarine}^{slope}`` | m``^2``/yr | 0.25 |
| Submarine slope diffusivity | ``\kappa_{submarine}^{slope}`` | m``^2``/yr | 100 |
| Submarine diffusion decay depth | ``\lambda_{submarine}`` | m | 1000 |
| Precipitation rate | ``R_{precip}`` | m/yr | 1.0 |
| Pelagic sedimentation rate | ``R_{pelagic}`` | mm/yr | 0.0 |
| Initial porosity of sediment at mudline | ``\phi_{o}`` | - | 0.4 |
| Sediment porosity decay depth | ``\lambda_{\phi}`` | m | 2500 |
#### tab:sediment-transport-example-parameters
*Example parameters for sediment transport model*