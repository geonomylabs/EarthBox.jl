using EarthBox
using Test

import EarthBox.Markers.MarkerTemperature.InitManager.HalfSpaceCooling: calculate_temperature
import EarthBox.ConversionFuncs: celsius_to_kelvin, cm_per_yr_to_meters_per_seconds

@testset "HalfSpaceCooling" begin
    surface_temperature_kelvin       = celsius_to_kelvin(0.0)
    asthenosphere_temperature_kelvin = celsius_to_kelvin(1300.0)
    distance_from_ridge_meters       = 1_000_000.0
    half_spreading_rate_m_s          = cm_per_yr_to_meters_per_seconds(5.0)
    thermal_diffusivity_m2_s         = 1.0e-6

    sample_depths_m = [0.0, 5_000.0, 25_000.0, 60_000.0, 100_000.0]

    expected_temperature_kelvin = [
        273.0,
        418.4922274291355,
        946.8557598815687,
        1454.356234340582,
        1566.6512337755787,
    ]

    for (k, depth_m) in enumerate(sample_depths_m)
        T = calculate_temperature(
            depth_m,
            surface_temperature_kelvin,
            asthenosphere_temperature_kelvin,
            distance_from_ridge_meters,
            half_spreading_rate_m_s,
            thermal_diffusivity_m2_s,
        )
        @test isapprox(T, expected_temperature_kelvin[k]; atol=1e-9)
    end
end
