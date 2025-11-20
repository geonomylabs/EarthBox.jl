module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .MaterialDescriptionGroup: MaterialDescription
import .RandomFrictionGroup: RandomFriction
import .SofteningGroup: Softening
import .ViscosityLimitsGroup: ViscosityLimits
import .StressLimitsYieldGroup: StressLimitsYield
import .StressLimitsPowerLawGroup: StressLimitsPowerLaw
import .HydrothermalGroup: Hydrothermal
import .CompactionGroup: Compaction
import .BoundaryFrictionGroup: BoundaryFriction
import .SerpentinizationGroup: Serpentinization
import .MeltDamageGroup: MeltDamage

"""
    Parameters <: AbstractParameterCollection

Parameter collection for material properties.

# Fields
- `material_description::`[`MaterialDescription`](@ref): Parameters for material model identifiers
- `random_friction::`[`RandomFriction`](@ref): Parameters for spatial and temporal variation of friction coefficient
- `softening::`[`Softening`](@ref): Parameters for strain softening behavior whereby pre-exponential is modified with increasing strain
- `viscosity_limits::`[`ViscosityLimits`](@ref): Minimum and maximum viscosity limits
- `stress_limits_yield::`[`StressLimitsYield`](@ref): Stress limits for plastic yielding
- `stress_limits_power_law::`[`StressLimitsPowerLaw`](@ref): Stress limits for power law creep
- `hydrothermal::`[`Hydrothermal`](@ref): Parameters for hydrothermal circulation model used to mimic the effects of hydrothermal circulation on thermal structure
- `compaction::`[`Compaction`](@ref): Parameters for sediment compaction
- `boundary_friction::`[`BoundaryFriction`](@ref): Parameters for boundary friction zones
- `serpentinization::`[`Serpentinization`](@ref): Parameters for serpentinization model
- `melt_damage::`[`MeltDamage`](@ref): Parameters for melt-induced damage

# Constructor
    Parameters()::Parameters

"""
mutable struct Parameters <: AbstractParameterCollection
    material_description::MaterialDescription
    random_friction::RandomFriction
    softening::Softening
    viscosity_limits::ViscosityLimits
    stress_limits_yield::StressLimitsYield
    stress_limits_power_law::StressLimitsPowerLaw
    hydrothermal::Hydrothermal
    compaction::Compaction
    boundary_friction::BoundaryFriction
    serpentinization::Serpentinization
    melt_damage::MeltDamage
end

function Parameters()::Parameters
    return Parameters(
        MaterialDescription(),
        RandomFriction(),
        Softening(),
        ViscosityLimits(),
        StressLimitsYield(),
        StressLimitsPowerLaw(),
        Hydrothermal(),
        Compaction(),
        BoundaryFriction(),
        Serpentinization(),
        MeltDamage()
    )
end

end # module 