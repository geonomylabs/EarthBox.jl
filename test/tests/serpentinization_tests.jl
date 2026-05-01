using EarthBox
using Test

import EarthBox.Serpentinization: calculate_temperature_dependent_rate_factor
import EarthBox.Serpentinization: calculate_plastic_strain_rate_rate_factor
import EarthBox.Serpentinization: calculate_serpentinization_rate
import EarthBox.Serpentinization: calculate_incremental_serpentinization_ratio
import EarthBox.ConversionFuncs: celsius_to_kelvin, years_to_seconds

@testset "Serpentinization" begin
    temperature_celsius                  = 400.0
    strain_rate_plastic                  = 1e-15
    timestep_seconds                     = years_to_seconds(50_000.0)
    nsteps                               = 3
    serpentinization_temperature_kelvins = celsius_to_kelvin(340.5)
    maximum_serpentinization_rate        = 1e-11
    nominal_strain_rate_serpentinization = 1e-13

    temperature_kelvins = celsius_to_kelvin(temperature_celsius)

    temperature_factor = calculate_temperature_dependent_rate_factor(
        temperature_kelvins, serpentinization_temperature_kelvins,
    )
    strain_rate_factor = calculate_plastic_strain_rate_rate_factor(
        strain_rate_plastic, nominal_strain_rate_serpentinization,
    )
    total_rate = calculate_serpentinization_rate(
        temperature_kelvins, strain_rate_plastic,
        maximum_serpentinization_rate,
        serpentinization_temperature_kelvins,
        nominal_strain_rate_serpentinization,
    )

    expected_temperature_factor = 0.009284732193196464
    expected_strain_rate_factor = 0.009950166250831893
    expected_total_rate         = 9.238462891675584e-16

    @test isapprox(temperature_factor, expected_temperature_factor; atol=1e-12)
    @test isapprox(strain_rate_factor, expected_strain_rate_factor; atol=1e-12)
    @test isapprox(total_rate,         expected_total_rate;         atol=1e-20)

    # Cumulative ratio after nsteps explicit-Euler increments.
    serpentinization_ratio = 0.0
    cumulative_after_step = zeros(nsteps)
    for istep in 1:nsteps
        increment = calculate_incremental_serpentinization_ratio(
            serpentinization_ratio, strain_rate_plastic, temperature_kelvins,
            timestep_seconds,
            serpentinization_temperature_kelvins,
            maximum_serpentinization_rate,
            nominal_strain_rate_serpentinization,
        )
        serpentinization_ratio += increment
        cumulative_after_step[istep] = serpentinization_ratio
    end

    expected_cumulative_after_step = [
        0.001457718582751707,
        0.0029133122220369144,
        0.0043667840154252005,
    ]
    for k in 1:nsteps
        @test isapprox(cumulative_after_step[k],
                       expected_cumulative_after_step[k]; atol=1e-12)
    end
end
