module CreepStrainRate

import ..ExponentialTerm: calculate_exponential_term

"""
    calculate_creep_strain_rate(
        pressure::Float64, temperature::Float64, gas_constant::Float64,
        pre_exponential::Float64, activation_energy::Float64, activation_volume::Float64,
        use_old_term::Bool, old_exponential_term::Float64, stress_invariant::Float64;
        stress_exponent::Float64 = 1.0
    )::Tuple{Float64, Float64}

Strain rate for dislocation or diffusion creep.

Dislocation or diffusion creep is modeled using the following formulation:

eii = A*exp_term*sii^n                                              (Eq. 1)

where A is the pre-exponential term, eii is the second invariant of the
strain rate tensor, sii is the second invariant of the stress tensor, n is
the power-law stress exponent (n = 1 for diffusion creep), and exp_term is
defined as:

exp_term = exp[-(Ea + Va*P)/RT)                                     (Eq. 2)

where A is the pre-exponential term, Ea is activation energy, Va is
activation volume, R is the gas constant and T is temperature.

# Arguments
- `pressure`: Marker pressure pressure (Pa)
- `temperature`: Marker temperature (K)
- `gas_constant`: Gas constant (J/mol/K)
- `pre_exponential`: Pre-exponential term (1/s/MPa^n)
- `activation_energy`: Activation energy (kJ/mol)
- `activation_volume`: Activation volume (J/MPa/mol) (cm^3)
- `stress_invariant`: Marker stress invariant in (Pa)
- `use_old_term`: Flag to use old terms. If this parameter is equal to true then the previously calculated exponential terms are reused
- `old_exponential_term`: Part of (Eq. 1) that does not depend on stress. This term can be reused from a previous calculation if use_old_term is true
- `stress_exponent`: Stress exponent n (for diffusion creep n = 1; use default value)

# Returns
- `strain_rate_invariant`: Second invariant of the strain rate tensor in 1/seconds
- `exponential_term`: Part of (Eq. 1) that does not depend on stress (see Eq. 2)
"""
@inline function calculate_creep_strain_rate(
    pressure_megapascals::Float64,
    temperature::Float64,
    gas_constant::Float64,
    pre_exponential::Float64,
    activation_energy::Float64,
    activation_volume::Float64,
    use_old_term::Bool,
    old_exponential_term::Float64,
    stress_invariant_megapascals::Float64;
    stress_exponent::Float64 = 1.0
)::Tuple{Float64, Float64}
    exponential_term = calculate_exponential_term(
        activation_energy, activation_volume, gas_constant,
        temperature, pressure_megapascals,
        use_old_term, old_exponential_term
    )
    # Try to avoid expensive power operations with floating point exponents
    if stress_exponent == 1.0
        stress_term = stress_invariant_megapascals
    else
        stress_term = stress_invariant_megapascals^stress_exponent
    end
    strain_rate_invariant = pre_exponential*exponential_term*stress_term
    return strain_rate_invariant, exponential_term
end

end # module 