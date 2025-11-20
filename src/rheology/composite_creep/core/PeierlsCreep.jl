module PeierlsCreep

using ..ExponentialTerm: calculate_activation_energy_term

const PaToMPa = 1.0/1e6
"""
    calculate_peierls_creep_strain_rate(
        pressure::Float64, temperature::Float64, gas_constant::Float64,
        pre_exponential::Float64, activation_energy::Float64,
        activation_volume::Float64, stress_exponent_m::Float64,
        stress_exponent_n::Float64, stress_peierls::Float64,
        use_old_term::Bool, old_activation_term::Float64,
        stress_invariant::Float64
    )::Tuple{Float64, Float64}

Calculate Peierls creep strain rate and activation term.

Peierls creep is modeled using the following formulation:

    eii = A*sii^2*exp[-activation_term*stress_inner_term]           (Eq. 1)

where eii is the second invariant of the strain rate tensor, sii is the
second invariant of the stress tensor, A is the pre-exponential term,
with the activation term defined as:

        activation_term = (Ea + Va*P)/RT                            (Eq. 2)

where Ea is activation energy, Va is activation volume, P is pressure,
R is the gas constant, T is temperature, and with the stress inner term
defined as:

    [1 - (sii/stress_peierls)^m]^n                                  (Eq. 3)

where stress_peierls is the stress for Peierls creep, m is the stress
exponent m, and n is the stress exponent n.

# Arguments
- `pressure`: Pressure in Pascals
- `temperature`: Temperature in Kelvin
- `gas_constant`: Gas constant in J/mol/K
- `pre_exponential`: Pre-exponential term in 1/s*1/MPa^2
- `activation_energy`: Activation energy in kJ/mol
- `activation_volume`: Activation volume in J/MPa/mol (cm^3)
- `stress_exponent_m`: Stress exponent m for non-dimensional exponential term
- `stress_exponent_n`: Stress exponent n for non-dimensional exponential term
- `stress_peierls`: Stress for Peierls creep in MPa
- `use_old_term`: Flag to use old terms (1 to use previously calculated activation term)
- `old_activation_term`: Activation term that does not depend on stress
- `stress_invariant`: Stress invariant in Pascals

# Returns
- `strain_rate_invariant`: Second invariant of the strain rate tensor in 1/s
- `activation_term`: Activation term that does not depend on stress
"""
@inline function calculate_peierls_creep_strain_rate(
    pressure_megapascals::Float64,
    temperature::Float64,
    gas_constant::Float64,
    pre_exponential::Float64,
    activation_energy::Float64,
    activation_volume::Float64,
    stress_exponent_m::Float64,
    stress_exponent_n::Float64,
    stress_peierls::Float64,
    use_old_term::Bool,
    old_activation_term::Float64,
    stress_invariant_megapascals::Float64
)::Tuple{Float64, Float64}
    activation_term = if use_old_term == false
        calculate_activation_energy_term(
            activation_energy, activation_volume, gas_constant,
            temperature, pressure_megapascals
        )
    else
        old_activation_term
    end
    stress_ratio = (stress_invariant_megapascals/stress_peierls)
    # Try to avoid expensive power operations with floating point exponents
    # stress_exponent_m is often 1 (avoid float exponents)
    if stress_exponent_m == 1.0
        stress_term_m = 1.0 - stress_ratio
    else
        stress_term_m = 1.0 - stress_ratio^stress_exponent_m
    end
    # stress_exponent_n is often 2 (avoid float exponents)
    if stress_exponent_n == 2.0
        stress_inner_term = stress_term_m*stress_term_m
    else
        stress_inner_term = stress_term_m^stress_exponent_n
    end
    exponential_term = exp(-activation_term*stress_inner_term)
    strain_rate_invariant = (
        pre_exponential
        * exponential_term
        * (stress_invariant_megapascals*stress_invariant_megapascals)
    )
    return strain_rate_invariant, activation_term
end

"""
    calculate_creep_stress_peierls(
        pressure::Float64,
        temperature::Float64,
        gas_constant::Float64,
        pre_exponential::Float64,
        activation_energy::Float64,
        activation_volume::Float64,
        use_old_term::Bool,
        old_activation_term::Float64,
        stress_peierls::Float64,
        strain_rate::Float64,
        stress_exponent_m::Float64,
        stress_exponent_n::Float64,
    )::Float64

Calculate the stress for Peierls creep for a given strain rate.

# Arguments
- `pressure`: Pressure in Pascals
- `temperature`: Temperature in Kelvin
- `gas_constant`: Gas constant in J/mol/K
- `pre_exponential`: Pre-exponential term in 1/s*1/MPa^2
- `activation_energy`: Activation energy in kJ/mol
- `activation_volume`: Activation volume in J/MPa/mol (cm^3)
- `iuse_old_term`: Flag to use old terms (1 to use previously calculated activation term)
- `old_activation_term`: Activation term that does not depend on stress
- `stress_peierls`: Stress for Peierls creep in MPa
- `strain_rate`: Strain rate in 1/s
- `stress_exponent_m`: Stress exponent m for non-dimensional exponential term
- `stress_exponent_n`: Stress exponent n for non-dimensional exponential term

# Returns
- `stress_invariant`: Stress invariant in Pascals
"""
function calculate_creep_stress_peierls(
    pressure::Float64,
    temperature::Float64,
    gas_constant::Float64,
    pre_exponential::Float64,
    activation_energy::Float64,
    activation_volume::Float64,
    use_old_term::Bool,
    old_activation_term::Float64,
    stress_peierls::Float64,
    strain_rate::Float64,
    stress_exponent_m::Float64,
    stress_exponent_n::Float64,
)::Float64
    nmax = 1000
    tolerance = 1e-5
    stress_left = 0.0  # Pa
    stress_right = stress_peierls*1e6  # Pa
    icount = 1
    criterion = 1e24
    stress_middle = 0.0

    while icount < nmax && criterion > tolerance
        strain_rate_left, _ = calculate_peierls_creep_strain_rate(
            pressure*PaToMPa,
            temperature,
            gas_constant,
            pre_exponential,
            activation_energy,
            activation_volume,
            stress_exponent_m,
            stress_exponent_n,
            stress_peierls,
            use_old_term,
            old_activation_term,
            stress_left*PaToMPa
        )
        delta_left = strain_rate - strain_rate_left
        stress_middle = stress_left + (stress_right - stress_left)/2.0
        strain_rate_middle, _ = calculate_peierls_creep_strain_rate(
            pressure*PaToMPa,
            temperature,
            gas_constant,
            pre_exponential,
            activation_energy,
            activation_volume,
            stress_exponent_m,
            stress_exponent_n,
            stress_peierls,
            use_old_term,
            old_activation_term,
            stress_middle*PaToMPa
        )
        delta_middle = strain_rate - strain_rate_middle
        criterion = (stress_right - stress_left)/2.0
        icount += 1
        if delta_left*delta_middle > 0.0
            stress_left = stress_middle
        else
            stress_right = stress_middle
        end
    end
    return stress_middle
end

end # module PeierlsCreep 