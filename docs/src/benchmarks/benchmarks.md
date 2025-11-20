# Benchmarks

Benchmarks can be run for each individual benchmark case or in batch mode.

For single execution and a description of each benchmark see 
[Single Benchmark Execution](@ref). For batch execution see 
[Batch Benchmark Execution](@ref).

!!! warning "execution time"
    Some benchmarks (e.g. couette flow and channel flow cases) may only take a few
    minutes to run on modern multithreaded hardware while other may take a half
    an hour to multiple hours.

The links provide more details about each benchmark packaged with EarthBox including
results:

- [Couette Flow with Viscous Heating](@ref)
- [Channel Flow Variable Conductivity](@ref)
- [Channel Flow Non-steady Temperature](@ref)
- [Channel Flow Non-Newtonian](@ref)
- [Solid Body Rotation](@ref)
- [Rayleigh Taylor Instability](@ref)
- [Elastic Slab](@ref)
- [Viscoelastic Stress Buildup](@ref)
- [Box Convection Isoviscous 1a](@ref)
- [Plasticity Benchmark Kaus10](@ref)
- [Seafloor Spreading](@ref)
- [Flexure Triangular Hole](@ref)
- [Viscoelastic Extension](@ref)
- [Viscoelastic Contraction](@ref)
- [Simple Sedimentation](@ref)

## Benchmark Model Inputs

All benchmarks packages with EarthBox are located in `src/benchmarks/models`. Inputs
are defined using yaml formatted input files `model.yaml` and `materials.yaml`.
The material library used for benchmarks is located at 
`src/material_library/registries/benchmarks/benchmarks.yml`.