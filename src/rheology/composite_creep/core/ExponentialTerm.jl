module ExponentialTerm

const KJ_TO_J = 1000.0

"""
    calculate_exponential_term(
        activation_energy::Float64,
        activation_volume::Float64,
        gas_constant::Float64,
        temperature::Float64,
        pressure_megapascals::Float64,
        iuse_old_exponential_term::Int64,
        exponential_term::Float64
    )::Float64

Calculate exponential term for dislocation creep or diffusion creep.

Exponential term is defined as:

    exponential_term = exp[-(Ea + Va*P)/RT]                               (Eq 1)

Ea is activation energy, Va is activation volume, P is pressure, R is the gas 
constant and T is temperature. The units of A are 1/s/MPa^n where n is the stress 
exponent (n = 1 for diffusion creep).

# Arguments
- `activation_energy`: Activation energy in kJ/mol.
- `activation_volume`: Activation volume in J/MPa/mol (cm^3).
- `gas_constant`: Gas constant in J/mol/K.
- `temperature`: Temperature in Kelvin.
- `pressure_megapascals`: Pressure in MPa.
- `iuse_old_exponential_term`: Flag to use old exponential terms. If this parameter is equal to 1 then
  the previously calculated exponential term is used.
- `exponential_term`: Part of the exponential term that does not depend on stress. This term
  can be reused from a previous calculation if iuse_old_term is 1.

# Returns
- `exponential_term`: Exponential term in 1/sMPa^n.
"""
@inline function calculate_exponential_term(
    activation_energy::Float64,
    activation_volume::Float64,
    gas_constant::Float64,
    temperature::Float64,
    pressure_megapascals::Float64,
    use_old_exponential_term::Bool,
    exponential_term::Float64
)::Float64
    if use_old_exponential_term == false
        activation_term = calculate_activation_energy_term(
            activation_energy, activation_volume, gas_constant,
            temperature, pressure_megapascals
        )
        exponential_term = exp(-activation_term)
    end
    return exponential_term
end

"""
    calculate_activation_energy_term(
        activation_energy::Float64,
        activation_volume::Float64,
        gas_constant::Float64,
        temperature::Float64,
        pressure_megapascals::Float64
    )::Float64

Calculate activation energy term for dislocation or diffusion creep.

Activation energy term is defined as:

    activation_term = (Ea + Va*P)/RT                                      (Eq 1)

where Ea is activation energy, Va is activation volume, P is pressure, R
is the gas constant and T is temperature.

# Arguments
- `activation_energy`: Activation energy in kJ/mol.
- `activation_volume`: Activation volume in J/MPa/mol (cm^3).
- `gas_constant`: Gas constant in J/mol/K.
- `temperature`: Temperature in Kelvin.
- `pressure_megapascals`: Pressure in MPa.

# Returns
- `activation_term`: Activation term (units cancel).
"""
@inline function calculate_activation_energy_term(
    activation_energy::Float64,
    activation_volume::Float64,
    gas_constant::Float64,
    temperature::Float64,
    pressure_megapascals::Float64
)::Float64
    fac1 = 1.0 / (gas_constant * temperature)
    activation_term = (
        activation_energy * KJ_TO_J + activation_volume * pressure_megapascals
        ) * fac1
    return activation_term
end

end # module ExponentialTerm 