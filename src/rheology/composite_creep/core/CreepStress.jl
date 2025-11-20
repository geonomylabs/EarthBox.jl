module Creep

import ..ExponentialTerm: calculate_exponential_term

"""
    calculate_creep_stress(
        pressure::Float64,
        temperature::Float64,
        gas_constant::Float64,
        pre_exponential::Float64,
        activation_energy::Float64,
        activation_volume::Float64,
        iuse_old_term::Int64,
        old_exponential_term::Float64,
        strain_rate::Float64;
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
- `strain_rate`: Strain rate invariant (1/s)
- `iuse_old_term`: Flag to use old terms. If this parameter is equal to 1 then the previously non-stress dependent terms are reused
- `old_exponential_term`: Part of (Eq. 1) that does not depend on stress. This term can be reused from a previous calculation if iuse_old_term is 1
- `stress_exponent`: Stress exponent n (for diffusion creep n = 1; use default value)

# Returns
- `stress_invariant`: Second invariant of the stress tensor Pa
- `exponential_term`: Part of (Eq. 1) that does not depend on stress (see Eq. 2)
"""
function calculate_creep_stress(
    pressure::Float64,
    temperature::Float64,
    gas_constant::Float64,
    pre_exponential::Float64,
    activation_energy::Float64,
    activation_volume::Float64,
    use_old_term::Bool,
    old_exponential_term::Float64,
    strain_rate::Float64;
    stress_exponent::Float64 = 1.0
)::Tuple{Float64, Float64}
    pressure_megapascals = pressure/1e6

    exponential_term = calculate_exponential_term(
        activation_energy, activation_volume, gas_constant,
        temperature, pressure_megapascals,
        use_old_term, old_exponential_term
    )

    stress_invariant_megapascals = (
        (
            strain_rate/pre_exponential/exponential_term
        )^(1.0/stress_exponent)
    )

    stress_invariant = stress_invariant_megapascals*1e6

    return stress_invariant, exponential_term
end

end # module 