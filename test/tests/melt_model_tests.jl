using EarthBox
using Random
using Test

import EarthBox.MeltModel.MeltFraction: get_melting_model_parameters
import EarthBox.MeltModel.MeltPropertiesOpt: calculate_delta_meltfrac_times_latent_heat_for_temperature
import EarthBox.MeltModel.MeltPropertiesOpt: calculate_delta_meltfrac_times_latent_heat_for_pressure
import EarthBox.MeltModel.MeltDamage: calculate_damage_factor_probabilistic
import EarthBox.MeltModel.MeltDamage: calculate_damage_probability
import EarthBox.MeltModel.MeltDamage: linear_melt_damage_probability_model
import EarthBox.MeltModel.Extraction: PartiallyMoltenZone
import EarthBox.FindShallowest: find_shallowest_partially_molten_mantle_marker_opt

const _MMT_DAMAGE_SEED = 20260430

@testset "MeltModel" begin

    @testset "PeridotiteSolidusLiquidus (Gerya2010)" begin
        itype_solidus = 3
        itype_liquidus = 3
        sample_pressures_pa = [
            0.0,
            3300 * 9.81 *  50_000.0,
            3300 * 9.81 * 150_000.0,
            3300 * 9.81 * 300_000.0,
            3300 * 9.81 * 600_000.0,
        ]
        expected_t_liquidus_k = [
            2073.0,
            2257.5261,
            2626.5783,
            3180.1566000000003,
            4287.3132000000005,
        ]
        expected_t_solidus_k = [
            1393.811,
            1595.5553443439599,
            1918.80830099564,
            2203.0984058825597,
            2502.8320922000003,
        ]
        for (k, p) in enumerate(sample_pressures_pa)
            t_liq, t_sol = get_melting_model_parameters(p, itype_solidus, itype_liquidus)
            @test isapprox(t_liq, expected_t_liquidus_k[k]; atol=1e-6)
            @test isapprox(t_sol, expected_t_solidus_k[k]; atol=1e-6)
        end
    end

    @testset "GabbroSolidusLiquidus (Gerya2010)" begin
        itype_solidus = 5
        itype_liquidus = 5
        sample_pressures_pa = [
            0.0,
            2800 * 9.81 *  2_000.0,
            2800 * 9.81 *  5_000.0,
            2800 * 9.81 *  8_000.0,
            2800 * 9.81 * 12_000.0,
        ]
        expected_t_liquidus_k = [
            1423.0,
            1428.76828,
            1437.4207,
            1446.07312,
            1457.60968,
        ]
        expected_t_solidus_k = [
            1327.0,
            1331.999176,
            1339.49794,
            1346.996704,
            1356.995056,
        ]
        for (k, p) in enumerate(sample_pressures_pa)
            t_liq, t_sol = get_melting_model_parameters(p, itype_solidus, itype_liquidus)
            @test isapprox(t_liq, expected_t_liquidus_k[k]; atol=1e-6)
            @test isapprox(t_sol, expected_t_solidus_k[k]; atol=1e-6)
        end
    end

    @testset "LatentHeat (Katz2003)" begin
        itype_solidus = 4
        itype_liquidus = 4
        latent_heat = 400.0e3      # J/kg
        density = 3300.0           # kg/m^3
        delta_temperature = 1.0    # K
        delta_pressure = 1000.0    # Pa
        T_geotherm_k = 1330.0 + 273.15

        sample_depths_m = [10_000.0, 50_000.0, 90_000.0, 120_000.0, 140_000.0]

        # Depths > 50 km are below the Katz2003 anhydrous solidus for this
        # geotherm, so the melt-fraction derivatives go to zero.
        expected_delta_heat_capacity = [
            600.6761265397072,
            715.566516643119,
            0.0,
            0.0,
            0.0,
        ]
        expected_delta_expansivity_x1e5 = [
            -12.806369806496432,
            -16.281742298270494,
            0.0,
            0.0,
            0.0,
        ]

        for (k, depth_m) in enumerate(sample_depths_m)
            pressure_pa = 3330 * 9.81 * depth_m
            dxLh_T = calculate_delta_meltfrac_times_latent_heat_for_temperature(
                pressure_pa, T_geotherm_k, delta_temperature,
                itype_solidus, itype_liquidus, latent_heat,
            )
            delta_heat_capacity = dxLh_T / (2.0 * delta_temperature)

            dxLh_P = calculate_delta_meltfrac_times_latent_heat_for_pressure(
                delta_pressure, pressure_pa, T_geotherm_k,
                itype_solidus, itype_liquidus, latent_heat,
            )
            delta_adiabatic = density * dxLh_P / (2.0 * delta_pressure)
            delta_expansivity_x1e5 = (delta_adiabatic / T_geotherm_k) * 1e5

            @test isapprox(delta_heat_capacity,
                           expected_delta_heat_capacity[k]; atol=1e-6)
            @test isapprox(delta_expansivity_x1e5,
                           expected_delta_expansivity_x1e5[k]; atol=1e-6)
        end
    end

    @testset "MeltDamage" begin
        avg_x = 250_000.0
        melt_damage_distance = 2_500.0
        melt_damage_factor = 10.0
        maximum_damage_probability = 0.8

        sample_x_offsets = [-2_500.0, -1_000.0, 0.0, 1_000.0, 2_500.0]
        expected_damage_probability = [
            0.0,
            0.523606797749979,
            0.8,
            0.523606797749979,
            0.0,
        ]
        for (k, dx) in enumerate(sample_x_offsets)
            x_marker = avg_x + dx
            p = calculate_damage_probability(
                x_marker, avg_x, melt_damage_distance, maximum_damage_probability,
            )
            @test isapprox(p, expected_damage_probability[k]; atol=1e-12)
        end

        # calculate_damage_factor_probabilistic uses Random.rand internally;
        # seed and pin the resulting sequence at the same x offsets above.
        Random.seed!(_MMT_DAMAGE_SEED)
        expected_damage_factor_probabilistic = [
            1.0,
            10.0,
            10.0,
            10.0,
            1.0,
        ]
        for (k, dx) in enumerate(sample_x_offsets)
            x_marker = avg_x + dx
            f = calculate_damage_factor_probabilistic(
                x_marker, avg_x, melt_damage_distance, melt_damage_factor,
                maximum_damage_probability,
            )
            @test isapprox(f, expected_damage_factor_probabilistic[k]; atol=1e-12)
        end

        # Linear damage-probability vs. crust height.
        h_threshold    = 500.0
        h_minimum      = 750.0
        h_intermediate = 2_000.0
        h_maximum      = 3_000.0
        intermediate_p = 0.1
        maximum_p      = 0.8
        sample_heights = [400.0, 700.0, 1_500.0, 2_500.0, 4_000.0]
        expected_linear_probability = [
            0.0,
            0.0,
            0.060000000000000005,
            0.45000000000000007,
            0.8,
        ]
        for (k, h) in enumerate(sample_heights)
            p = linear_melt_damage_probability_model(
                h, h_threshold, h_minimum, h_intermediate, h_maximum,
                maximum_p, intermediate_p,
            )
            @test isapprox(p, expected_linear_probability[k]; atol=1e-12)
        end
    end

    @testset "MagmaBodyExtraction" begin
        nmarkers_x = 101
        nmarkers_y = 101
        dx_m = 1.0
        dy_m = 1.0
        matid_background = Int16(1)
        matid_pm         = Int16(2)
        matid_magma      = Int16(3)
        magma_fraction   = 0.3
        mantle_melting_mat_ids = Int16[Int16(2)]

        nmarkers = nmarkers_x * nmarkers_y
        marker_x = zeros(Float64, nmarkers)
        marker_y = zeros(Float64, nmarkers)
        marker_matid = fill(matid_background, nmarkers)
        idx = 1
        for j in 1:nmarkers_y, i in 1:nmarkers_x
            marker_x[idx] = (i - 1) * dx_m
            marker_y[idx] = (j - 1) * dy_m
            idx += 1
        end

        # Box that flags partial-melt markers (matches demo fractions).
        x_max = maximum(marker_x); y_max = maximum(marker_y)
        xpm_o = x_max * 0.25; xpm_f = x_max * 0.75
        ypm_o = y_max * 0.10; ypm_f = y_max * 0.35
        partial_melt_flags_tmp = zeros(Int64, nmarkers)
        npm = 0
        for i in 1:nmarkers
            if xpm_o <= marker_x[i] <= xpm_f && ypm_o <= marker_y[i] <= ypm_f
                marker_matid[i] = matid_pm
                npm += 1
                partial_melt_flags_tmp[npm] = i
            end
        end
        partial_melt_flags = partial_melt_flags_tmp[1:npm]

        nmarkers_magma = floor(Int, magma_fraction * npm)

        nlayers = 10
        layer_counts = zeros(Int, nlayers)
        layer_offsets = zeros(Int, nlayers + 1)
        layered_partial_melt_indices = fill(Int64(-1), nmarkers)

        PartiallyMoltenZone.construct_layered_partially_molten_arrays!(
            marker_x, marker_y, marker_matid, npm,
            mantle_melting_mat_ids, partial_melt_flags,
            layer_counts, layer_offsets, layered_partial_melt_indices,
        )

        for _ in 1:nmarkers_magma
            (imarker_shallow, _y) = find_shallowest_partially_molten_mantle_marker_opt(
                marker_matid, marker_y, mantle_melting_mat_ids,
                layer_counts, layer_offsets, layered_partial_melt_indices,
            )
            if imarker_shallow != -999
                marker_matid[imarker_shallow] = matid_magma
            end
        end

        expected_npm                   = 1326
        expected_nmarkers_magma_target = 397
        expected_actual_magma_count    = 397
        expected_actual_pm_remaining   = 929
        expected_matid_checksum        = 11924

        @test npm == expected_npm
        @test nmarkers_magma == expected_nmarkers_magma_target
        @test count(==(matid_magma), marker_matid) == expected_actual_magma_count
        @test count(==(matid_pm),    marker_matid) == expected_actual_pm_remaining
        @test sum(Int.(marker_matid)) == expected_matid_checksum
    end

end
