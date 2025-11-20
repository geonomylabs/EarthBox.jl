module CompositeCreep

include("core/ExponentialTerm.jl")
include("core/CreepStrainRate.jl")
include("core/CreepStress.jl")
include("core/PeierlsCreep.jl")

using .CreepStrainRate: calculate_creep_strain_rate
using .PeierlsCreep: calculate_peierls_creep_strain_rate

""" Calculate composite viscosity and strain rate including:
- Dislocation Creep (Power Law Creep)
- Peierls Creep (Exponential Creep)
- Diffusion Creep

This function allows for re-use of terms (i.e. old terms) that do not depend on stress.
"""
@inline function calculate_composite_creep_viscosity_and_strain_rate(
        pressure_megapascals::Float64,
        temperature_kelvin::Float64,
        gas_constant::Float64,
        pre_exponential_dislocation::Float64,
        activation_energy_dislocation::Float64,
        activation_volume_dislocation::Float64,
        stress_exponent_dislocation::Float64,
        pre_exponential_peierls::Float64,
        stress_exponent_peierls_m::Float64,
        stress_exponent_peierls_n::Float64,
        stress_peierls::Float64,
        pre_exponential_diffusion::Float64,
        activation_energy_diffusion::Float64,
        activation_volume_diffusion::Float64,
        stress_invariant_megapascals::Float64,
        use_old_term::Bool,
        exponential_term_dislocation::Float64,
        activation_energy_term_peierls::Float64,
        exponential_term_diffusion::Float64
)::Tuple{Float64, Float64, Float64, Float64, Float64}

    # Dislocation Creep
    (
        strain_rate_invariant_dislocation,
        exponential_term_dislocation
    ) = calculate_creep_strain_rate(
        pressure_megapascals,
        temperature_kelvin,
        gas_constant,
        pre_exponential_dislocation,
        activation_energy_dislocation,
        activation_volume_dislocation,
        use_old_term,
        exponential_term_dislocation,
        stress_invariant_megapascals;
        stress_exponent = stress_exponent_dislocation
    )

    strain_rate_invariant_peierls = 0.0
    if use_peierls_creep(pre_exponential_peierls, temperature_kelvin)
        (
            strain_rate_invariant_peierls,
            activation_energy_term_peierls
        ) = calculate_peierls_creep_strain_rate(
            pressure_megapascals,
            temperature_kelvin,
            gas_constant,
            pre_exponential_peierls,
            activation_energy_dislocation,
            activation_volume_dislocation,
            stress_exponent_peierls_m,
            stress_exponent_peierls_n,
            stress_peierls,
            use_old_term,
            activation_energy_term_peierls,
            stress_invariant_megapascals
        )
    end

    strain_rate_invariant_diffusion = 0.0
    if use_diffusion_creep(pre_exponential_diffusion)
        (
            strain_rate_invariant_diffusion,
            exponential_term_diffusion
        ) = calculate_creep_strain_rate(
            pressure_megapascals,
            temperature_kelvin,
            gas_constant,
            pre_exponential_diffusion,
            activation_energy_diffusion,
            activation_volume_diffusion,
            use_old_term,
            exponential_term_diffusion,
            stress_invariant_megapascals;
            stress_exponent = 1.0
        )
    end

    strain_rate_invariant = calculate_effective_strain_rate_invariant(
        strain_rate_invariant_dislocation,
        strain_rate_invariant_peierls,
        strain_rate_invariant_diffusion
    )

    viscosity = calculate_effective_viscosity(
        stress_invariant_megapascals*1e6,
        strain_rate_invariant
    )

    return (
        viscosity,
        strain_rate_invariant,
        exponential_term_dislocation,
        activation_energy_term_peierls,
        exponential_term_diffusion
    )
end

@inline function use_diffusion_creep(pre_exponential_diffusion::Float64)::Bool
    return pre_exponential_diffusion > 0.0
end

""" Check if Peierls creep should be used.

Peierls creep is only used if the pre-exponential term is greater than zero
and the temperature is less than 1200 K.
"""
@inline function use_peierls_creep(pre_exponential_peierls::Float64, temperature::Float64)::Bool
    return pre_exponential_peierls > 0.0 && temperature < 1200.0
end

@inline function calculate_effective_strain_rate_invariant(
        strain_rate_invariant_dislocation::Float64,
        strain_rate_invariant_peierls::Float64,
        strain_rate_invariant_diffusion::Float64
)::Float64
    return max(strain_rate_invariant_dislocation, strain_rate_invariant_peierls) + 
           strain_rate_invariant_diffusion
end

@inline function calculate_effective_viscosity(
        stress_invariant_pascals::Float64,
        strain_rate_invariant::Float64
)::Float64
    return 0.5 * stress_invariant_pascals / strain_rate_invariant
end

end # module CompositeCreep 