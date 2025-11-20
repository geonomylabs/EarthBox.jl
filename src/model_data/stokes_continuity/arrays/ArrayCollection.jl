module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import EarthBox.Arrays.ArrayTypes.RhsStokesArray1D: RhsStokesArray1DState
import .BasicGridVelocityGroup: BasicGridVelocity
import .DensityGroup: Density
import .PlasticDefGroup: PlasticDef
import .PressureGroup: Pressure
import .ResidualStokesGroup: ResidualStokes
import .RhsStokesGroup: RhsStokes
import .ShearModulusGroup: ShearModulus
import .StaggeredGridVelocityGroup: StaggeredGridVelocity
import .StokesSolutionGroup: StokesSolution
import .StrainRateAndSpinGroup: StrainRateAndSpin
import .StressGroup: Stress
import .StressChangeGroup: StressChange
import .ViscosityGroup: Viscosity
import .VelocitySolutionGroup: VelocitySolution

"""
    Arrays <: AbstractArrayCollection

Data structure containing array groups for Stokes-continuity solver.

# Fields
- `rhs::`[`RhsStokes`](@ref): Right-hand side arrays for Stokes equations
- `staggered_grid_velocity::`[`StaggeredGridVelocity`](@ref): Velocity on staggered grid
- `basic_grid_velocity::`[`BasicGridVelocity`](@ref): Velocity on basic grid
- `velocity_solution::`[`VelocitySolution`](@ref): Velocity solution arrays
- `stokes_solution::`[`StokesSolution`](@ref): Combined Stokes solution arrays
- `strain_rate_and_spin::`[`StrainRateAndSpin`](@ref): Strain rate and spin tensor
- `stress::`[`Stress`](@ref): Stress tensor components
- `stress_change::`[`StressChange`](@ref): Stress change arrays
- `pressure::`[`Pressure`](@ref): Pressure field arrays
- `density::`[`Density`](@ref): Density field arrays
- `shear_modulus::`[`ShearModulus`](@ref): Shear modulus arrays
- `viscosity::`[`Viscosity`](@ref): Viscosity field arrays
- `plastic_def::`[`PlasticDef`](@ref): Plastic deformation arrays
- `residuals::`[`ResidualStokes`](@ref): Residual arrays for Stokes equations

# Constructor
    Arrays(ynum::Int, xnum::Int)::Arrays

Create a new Arrays collection with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

"""
mutable struct Arrays <: AbstractArrayCollection
    rhs::RhsStokes
    staggered_grid_velocity::StaggeredGridVelocity
    basic_grid_velocity::BasicGridVelocity
    velocity_solution::VelocitySolution
    stokes_solution::StokesSolution
    strain_rate_and_spin::StrainRateAndSpin
    stress::Stress
    stress_change::StressChange
    pressure::Pressure
    density::Density
    shear_modulus::ShearModulus
    viscosity::Viscosity
    plastic_def::PlasticDef
    residuals::ResidualStokes
end

function Arrays(ynum::Int, xnum::Int)::Arrays
    return Arrays(
        RhsStokes(ynum, xnum),
        StaggeredGridVelocity(ynum, xnum),
        BasicGridVelocity(ynum, xnum),
        VelocitySolution(ynum, xnum),
        StokesSolution(ynum, xnum),
        StrainRateAndSpin(ynum, xnum),
        Stress(ynum, xnum),
        StressChange(ynum, xnum),
        Pressure(ynum, xnum),
        Density(ynum, xnum),
        ShearModulus(ynum, xnum),
        Viscosity(ynum, xnum),
        PlasticDef(ynum, xnum),
        ResidualStokes(ynum, xnum)
    )
end

end # module 