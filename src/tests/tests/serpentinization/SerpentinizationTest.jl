module SerpentinizationTest

import EarthBox.Serpentinization: calculate_temperature_dependent_rate_factor
import EarthBox.Serpentinization: calculate_plastic_strain_rate_rate_factor
import EarthBox.Serpentinization: calculate_serpentinization_rate
import EarthBox.Serpentinization: calculate_incremental_serpentinization_ratio
import EarthBox.ConversionFuncs: celsius_to_kelvin, years_to_seconds

function run_test()::Nothing
    temperature_celsius = 400.0
    strain_rate_plastic = 1e-15
    timestep_years = 50_000.0
    nsteps = 3
    serpentinization_temperature_kelvins = celsius_to_kelvin(340.5)
    maximum_serpentinization_rate = 1e-11
    nominal_strain_rate_serpentinization = 1e-13

    time_step_seconds = years_to_seconds(timestep_years)
    temperature_kelvins = celsius_to_kelvin(temperature_celsius)

    temperature_factor = calculate_temperature_dependent_rate_factor(
        temperature_kelvins, serpentinization_temperature_kelvins
    )
    println("temp_factor: ", temperature_factor)

    strain_rate_factor = calculate_plastic_strain_rate_rate_factor(
        strain_rate_plastic,
        nominal_strain_rate_serpentinization
    )
    println("strain_rate_factor: ", strain_rate_factor)

    total_rate = calculate_serpentinization_rate(
        temperature_kelvins,
        strain_rate_plastic,
        maximum_serpentinization_rate,
        serpentinization_temperature_kelvins,
        nominal_strain_rate_serpentinization
    )
    println("total_rate: ", total_rate)

    serpentinization_ratio = 0.0
    for istep in 1:nsteps
        incremental_serpentinization_ratio = calculate_incremental_serpentinization_ratio(
            serpentinization_ratio, 
            strain_rate_plastic,
            temperature_kelvins, 
            time_step_seconds,
            serpentinization_temperature_kelvins,
            maximum_serpentinization_rate,
            nominal_strain_rate_serpentinization
        )
        serpentinization_ratio += incremental_serpentinization_ratio
        println("istep, serpentinization_ratio: ", istep, " ", serpentinization_ratio)
    end

    return nothing
end

end # module 