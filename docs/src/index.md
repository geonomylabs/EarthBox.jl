# Introduction to EarthBox

EarthBox is a multiphase visco-elasto-plastic marker-in-cell geodynamic modeling program that 
discretizes the Stokes-continuity and heat transport equations on a staggered grid using 
conservative finite differences and free-surface stabilization (See [Methods](@ref)). 
Advective processes are modeled using a 4-th order Runge-Kutta scheme that includes both 
grid and sub-grid changes in temperature and deviatoric stress. A variety of geologic 
processes are integrated in the EarthBox modeling system including partial melting, 
instantaneous melt extraction and transport, gabbroic melt fractionation, melt extrusion, 
lava flow, submarine hydrothermal cooling, exothermic serpentinization and sediment 
transport with erosion, deposition, compaction and an isostatically realistic global sea level.
EarthBox uses a composite viscous flow law that includes dislocation creep, diffusion creep,
Peierls creep, the effects of viscous softening associated with grain size reduction and the 
effects of partial melting. The frictional plastic model used by EarthBox includes strain weakening
and melt damage associated with melt transport networks above zones of melt focusing.

The goals of the EarthBox project are two-fold: (1) create an easy-to-use and easy-to-modify 
geodynamic modeling tool tuned for shared memory machines and small clusters with performance 
comparable to C++ and Fortran and (2) develop a system that can be used to build
easily extensible libraries of options and model cases. EarthBox is currently 2D with an experimental 
3D multigrid solver that will be integrated into the EarthBox system during a future release.
EarthBox is developed using the Julia programming language and includes arm's-length integration with the 
parallel direct solver MUMPS [amestoy01, amestoy19](@cite) that manages the stochastic 
nature of this powerful parallel direct solver that if not properly addressed can
lead to failed simulation runs.

To get started with EarthBox see [EarthBox Installation](@ref), 
[Building and Executing Models with the EarthBox API](@ref), [EarthBox Models](@ref) and
[EarthBox API Overview](@ref). 

The fundamental aspects of the marker-in-cell approach for solving incompressible flows with a 
free-surface used by EarthBox were first published by [harlow65](@citet). This approach was later 
modified for visco-elasto-plastic problems by [gerya03, gerya07, gerya2010](@citet). EarthBox uses a 
visco-elasto-plastic thermo-mechanical algorithm that closely follows the algorithm described 
by [gerya2010](@citet) but with the node-based plasticity approach from [gerya2019](@citet) to 
improve the convergence of the non-linear visco-elasto-plastic Stokes-continuity solver 
(see [Methods](@ref) and [Main Steps](@ref)).

EarthBox includes a melting model with linearized melt fraction and instantaneous melt 
extraction with optional emplacement at the Moho to simulate magma intrusion at mid-ocean
ridges (see [Melt generation, Transport and Emplacement](@ref)). Novel aspects of EarthBox's melting 
model are (1) a simple gabbroic fractionation model based on the proximity of extracted melt 
to the Moho, (2) a melt extrusion and lava flow model based on cellular automata that 
allows the modeling of the formation of seaward dipping reflectors
(see [Melt Extrusion and Lava Flow](@ref)) and (3) a probabilistic melt-damage model that 
weakens frictional plastic properties above regions of melt focusing that plays a critical role
in stabilizing spreading centers when extrusive processes are included.

Topography can be difficult to accurately track with particle-in-cell models. EarthBox
implements an Eulerian-Lagrangian approach whereby a topography marker chain is
advected using a 4th order Runge-Kutta scheme and then topographic elevation is 
interpolated to an Eulerian topography grid. The approach led to highly accurate tracking of
the air-rock interface comparable to finite-element Eulerian-Lagrangian codes. 

EarthBox provides an option to define global sea level using the isostatic equilibrium 
between a reference lithospheric column and the average pressure at the base of the 
model [kneller23](@cite). Sediment transport including erosion and deposition with fluvial 
and marine processes is implemented using the approach of [martinez19](@citet) 
(see [Sediment Transport](@ref)) but with the added effect of sediment compaction. 

EarthBox also includes the effects of hydrothermal circulation on submarine thermal 
structure that takes into account the effects of sediment overburden 
(see [Hydrothermal Circulation](@ref)) and how exothermic serpentinization reactions 
impact thermal conductivity and density (see [Serpentinization](@ref)).

Numerical and analytical benchmarks of EarthBox are presented in the [Benchmarks](@ref) providing a 
direct comparison to similar codes such as I2ELVIS that has been used 
extensively in published geodynamic modeling studies [gerya03, gerya2010, gerya2019](@cite).
EarthBox includes an integrated benchmarking system allowing the user to quickly
reproduce benchmarks, gain confidence in multiple aspects of the code and test the code 
after modifications are made (see [Benchmarks](@ref)).
Some novel aspects of the design of the EarthBox software are as follows:

- A consistent option scheme utilizing a common option state and Julia's multiple dispatch system
   making it easy to customize and extend EarthBox and manage complex and diverse model options.
- Integrated boundary-condition options for the Stokes-continuity and heat equations
   for specific problem types that take marker recycling into account (see [Boundary Condition Models](@ref)).
- A self-contained fully integrated benchmarking system allows most aspects of the code to be
   tested in a fully integrated manner (see [Benchmarks](@ref)).
- Flexible model building, execution and plotting features allowing the use of an API
   or YAML input files (see [EarthBox API Overview](@ref)).
- A multiple model case management system (see [Case Management](@ref)).





