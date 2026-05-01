using EarthBox
using Test

import EarthBox.RockProperties.RhoCpModel.UpdateManager.TemperatureDependentWaples:
    calculate_heat_capacity_waples
import EarthBox.RockProperties.ThermalConductivityModel.UpdateManager.Liao14:
    thermal_conductivity_liao14
import EarthBox.RockProperties.ThermalConductivityModel.UpdateManager.SekiguchiWaples:
    thermal_conductivity_sekiguchi_waples
import EarthBox.ConversionFuncs: celsius_to_kelvin

@testset "RockProperties" begin

    @testset "HeatCapacity (Waples T-dependent)" begin
        heat_capacity_at_20celsius = 800.0
        sample_temperatures_celsius = [0.0, 200.0, 500.0, 1000.0, 1330.0]

        expected_cp = [
            762.4,
            1045.7024000000001,
            1230.5,
            1279.1999999999998,
            1428.6873335999999,
        ]
        for (k, T_c) in enumerate(sample_temperatures_celsius)
            cp = calculate_heat_capacity_waples(T_c, heat_capacity_at_20celsius)
            @test isapprox(cp, expected_cp[k]; atol=1e-9)
        end
    end

    @testset "ThermalConductivity (Liao14)" begin
        # Ultra-basic rocks from Clauser and Huenges (1995) Table 1.
        k0 = 0.73
        a  = 1293.0
        density = 3300.0
        gravity = 9.8

        sample_depths_m = [0.0, 10_000.0, 30_000.0, 60_000.0, 100_000.0]
        # Linear geotherm 0..1330 °C over 0..160 km, matching the demo.
        temperature_gradient_C_per_km = 1330.0 / 160.0

        expected_k_liao = [
            4.424285714285714,
            3.7636544678180903,
            3.0014980043598825,
            2.4352846215774564,
            2.0765804466820126,
        ]
        for (idx, depth_m) in enumerate(sample_depths_m)
            T_k = celsius_to_kelvin(0.0) +
                  temperature_gradient_C_per_km * (depth_m / 1000.0)
            P_pa = density * gravity * depth_m
            k_val = thermal_conductivity_liao14(k0, a, P_pa, T_k)
            @test isapprox(k_val, expected_k_liao[idx]; atol=1e-9)
        end
    end

    @testset "ThermalConductivity (Sekiguchi-Waples)" begin
        # Peridotite (typical) from Hantschel and Kauerauf (2009) Table E.4.
        k0 = 4.0
        sample_temperatures_celsius = [0.0, 200.0, 500.0, 1000.0, 1330.0]

        expected_k_sw = [
            4.198811339721612,
            2.9740664721860464,
            2.3252527292936613,
            1.923460505690495,
            1.7955838016119776,
        ]
        for (k, T_c) in enumerate(sample_temperatures_celsius)
            T_k = celsius_to_kelvin(T_c)
            k_val = thermal_conductivity_sekiguchi_waples(k0, T_k)
            @test isapprox(k_val, expected_k_sw[k]; atol=1e-9)
        end
    end
end
