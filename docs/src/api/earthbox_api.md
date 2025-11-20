# EarthBox API Overview

The EarthBox API provides access to model initialization and execution functions.

See [Building and Executing Models with the EarthBox API](@ref) for examples
of how to use the EarthBox API for building models and running time steps.

## API Components

- [Initializing EarthBox](@ref)
    - Earthbox models are initialized by creating a `EarthBoxState` object
- [Staggered Grid Initialization](@ref)
    - A staggered grid can be initialized using a variety of grid types and
       refinement parameters.
- A variety of option are available for defining material geometries including
   [StickyAirGeometry](@ref), , [EarthLayering](@ref), [LithoStrongZones](@ref), 
   [Plume](@ref), [CrustalHole](@ref), [RayleighTaylor](@ref),
   [MobileWall](@ref), [Sandbox](@ref), [InternalVelocityZone](@ref), 
   [FractureZone](@ref), [WeakFault](@ref) and [WeakSeed](@ref)
- [BoundaryConditions](@ref)
    - A variety of integrated boundary condition models are available for the
       Stokes-continuity and heat equations that include marker recycling
       schemes. Initialization API functions are available for pressure, 
       temperature, velocity, velocity stepping, velocity stopping, defining
       velocity with strain rate and for defining transient bottom temperature 
       boundary conditions that approximate the thermal effects of a plume 
       spreading out beneath a tectonic plate.
- [MarkerCoordinates](@ref)
    - Markers can be initialized using a regular or random distribution.
- [MarkerTemperature](@ref)
    - Marker temperature can be initialized using a variety of models including
       analytical models for the Earth's geotherm.
- [Marker Friction Coefficients](@ref)
    - Friction coefficients defined in markers can be initialized using a
       a regular or randomized model, and time-dependent randomization can
       be activated.
- [MarkerCohesion](@ref)
    - Initialize marker cohesion. 
- [MarkerDilatancy](@ref)
    - Initialize marker dilatancy.
- [MarkerViscousStrainSoftening](@ref)
    - A viscous strain softening model can be implemented whereby the 
       pre-exponential terms of dislocation creep are modified by a function
       of plastic strain.
- [MarkerStressLimits](@ref)
    - Initialize stress limits for power-law rheology and plastic yield stress.
- [MarkerBoundaryFriction](@ref)
    - Initialize boundary friction layer thickness, friction angle and cohesion.
- [MarkerMaterials](@ref)
    - Initialize a geometric marker material initialization model, define
       viscosity limits, and activate the use of fluid pressure in yield stress
       equations, and define a plastic healing rate.
- [RockProperties](@ref)
    - Initialize rock property models for density, thermal conductivity and
       heat capacity.
- [Hydrothermal Circulation](@ref)
    - Initialize a hydrothermal circulation model that approximates the effects 
       of hydrothermal circulation on thermal structure at spreading centers using
       effective thermal conductivity.
- [Serpentinization Model](@ref)
    - Initialize a serpentinization model that affects thermal conductivity and 
       density and produces heat during exothermic serpentinization reactions.
- [StokesContinuitySolver](@ref)
    - Initialize the Stokes-continuity including gravitational acceleration, 
       activation of interface stabilization to avoid the drunken sailor effect,
       an option for turning off gravity at a specified time and an option to 
       turn off the solver or use kinematic models instead of solving for velocity.
- [GlobalPlasticityLoop](@ref)
    - Initialize the global plasticity loop (Picard iterations) used to solve 
       the non-linear visco-elasto-plastic Stoke-continuity equations including 
       tolerance, maximum number of iterations and the type of global plasticity
       loop with option for defining plastic failure on nodes or markers.
- [HeatSolver](@ref)
    - Initialize the heat-equation solver including options for shear heating, 
       adiabatic heating, maximum temperature change used for adaptive time stepping,
       and a sticky air thermal correction that increases the accuracy of thermal 
       solutions at the sticky-rock interface.
- [Marker Advection](@ref "Advection")
    - Initialize the marker advection scheme including adaptive time stepping
       parameters and sub-grid thermal and stress diffusion.
- [Interpolation](@ref)
    - Initialization of options involving interpolation including harmonic
       averaging of shear viscosity on basic nodes to define viscosity on
       pressure nodes.
- [MeltModel](@ref)
    - Initialize the melting model include melt fraction calculation using
       specified solidus and liquidus models and melt-dependent rock property 
       models (viscosity, density and heat capacity).
- [MeltExtraction](@ref)
    - Initialize the melt extraction model for extracting partial melt from
        fertile ultramafic rocks with option for injecting magma at the base of
        the moho to model rapid melt transport and a gabbro fractionation
        model that transforms the properties of gabbroic melt based on
        proximity to the Moho.       
- [MeltExtrusion](@ref)
    - Initialize a melt extrusion model that models lava flow using a cellular 
       automata approach.
- [MeltDamage](@ref)
    - Initialize a probabilistic melt damage model hat modifies friction angles
       to approximate the weakening effect of melt transport networks above
       melt focusing points in the partially molten mantle.
- [Topography](@ref)
    - Initialize a topography model that uses a marker-chain with a
       Lagrangian-Eulerian approach for tracking the sticky-rock interface.
- [Sealevel](@ref)
    - Initialize a sea level model with different models for defining sea level
       including constant, using the upper-left corning the the lithosphere and
       using isostatic balance between a reference lithosphere column and average 
       pressure at the base of the model domain. A base level shift model can 
       also be defined.
- [RelativeBaseLevel](@ref)
    - Initialize the reference lithosphere used in the seal-level model based
       on average pressure at the base of the model. 
- [SedimentTransport](@ref)
    - Initialize the sediment transport model with deposition, erosion, fluvial
       transport and sediment compaction.
- [Salt Deposition](@ref)
    - Initialize a salt deposition model.
- [Time Step Execution](@ref)
    - Initialize time stepping and velocity stepping parameters and execute
       model times steps.
- [Case Management](@ref)
    - Define a dictionary of parameters for multiple cases and define model
       parameters based on a case name.






