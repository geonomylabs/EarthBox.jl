using EarthBox
using Test

import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel.DensityProps:
    DensityModel, DensityProperties
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel:
    calculate_isostatic_relative_base_level
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel.ReferenceLithosphere:
    LithosphereThicknesses, LithosphereThermalProps, reference_lithosphere
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel.LithostaticPressure:
    calculate_lithostatic_pressure_from_marker_swarm
import EarthBox.MathTools: linear_interp

# Density model with expansivity = compressibility = 0 (so the analytical
# isostatic formula is exact). Crust = 2800, mantle = 3300.
function _bl_simple_density_model()::DensityModel
    upr_crust = DensityProperties(2800.0, 0.0, 0.0)
    lwr_crust = DensityProperties(2800.0, 0.0, 0.0)
    upr_lith  = DensityProperties(3300.0, 0.0, 0.0)
    mid_lith  = DensityProperties(3300.0, 0.0, 0.0)
    lwr_lith  = DensityProperties(3300.0, 0.0, 0.0)
    asth      = DensityProperties(3300.0, 0.0, 0.0)
    return DensityModel(upr_crust, lwr_crust, upr_lith, mid_lith, lwr_lith, asth)
end

function _bl_marker_swarm(xnum, ynum, mpcx, mpcy, xsize, ysize)
    nx = xnum * mpcx
    ny = ynum * mpcy
    nmarkers = nx * ny
    marker_x = zeros(nmarkers)
    marker_y = zeros(nmarkers)
    dx = xsize / xnum / mpcx
    dy = ysize / ynum / mpcy
    k = 1
    for j in 1:ny, i in 1:nx
        marker_x[k] = (i - 1) * dx
        marker_y[k] = (j - 1) * dy
        k += 1
    end
    return marker_x, marker_y, dx, dy
end

function _bl_marker_density(gridy, density_gridy, marker_y, sticky_thickness)
    nnodes = length(gridy)
    nmarkers = length(marker_y)
    marker_rho = zeros(nmarkers)
    for i in 1:nmarkers
        y = marker_y[i] - sticky_thickness
        marker_rho[i] = linear_interp(
            y, gridy, density_gridy, nnodes, gridy[1], gridy[end], nnodes,
        )
    end
    return marker_rho
end

@testset "RelativeBaseLevel" begin
    xsize = 140_000.0
    ysize = 140_000.0
    xnum  = 160
    ynum  = 160
    sticky_thickness = 10_000.0
    mpcx = 8
    mpcy = 8

    density_model = _bl_simple_density_model()

    iuse_linear_segments = 0
    model_lith_thickness = LithosphereThicknesses(
        ysize - sticky_thickness, 22_000.0, 10_000.0, 122_000.0, 100.0,
    )
    model_lith_thermal = LithosphereThermalProps(
        0.0, 600.0, 1330.0, 0.4, 2.5, 2.5, 2.25, 9e-7, 9e-7, 0.0, 125_000.0,
    )

    # Build the model lithosphere temperature/density/pressure profile.
    gridy, temp_gridy, density_gridy, pressure_gridy = reference_lithosphere(
        density_model, model_lith_thickness, model_lith_thermal;
        number_of_pressure_iterations=3,
        iuse_linear_segments=iuse_linear_segments,
    )

    # Marker swarm + interpolated marker densities.
    marker_x, marker_y, marker_dx, marker_dy = _bl_marker_swarm(
        xnum, ynum, mpcx, mpcy, xsize, ysize)
    marker_rho = _bl_marker_density(gridy, density_gridy, marker_y, sticky_thickness)

    # Production lithostatic-pressure-from-markers reduction.
    nmarkers = length(marker_x)
    marker_x_scratch   = zeros(nmarkers)
    marker_y_scratch   = zeros(nmarkers)
    marker_rho_scratch = zeros(nmarkers)
    y_topo = sticky_thickness
    (
        gridy_marker_avg, density_gridy_from_markers, pressure_gridy_from_markers,
    ) = calculate_lithostatic_pressure_from_marker_swarm(
        marker_x_scratch, marker_y_scratch, marker_rho_scratch,
        marker_x, marker_y, marker_rho,
        ysize, marker_dx * 8, y_topo, marker_dy * 4, marker_dx * 4,
    )
    pressure_at_base_from_markers = pressure_gridy_from_markers[end]

    # Reference lithosphere (matching crust thickness, fuller column).
    reference_lith_thickness = LithosphereThicknesses(
        ysize, 22_000.0, 10_000.0, 122_000.0, 10.0,
    )
    reference_lith_thermal = LithosphereThermalProps(
        0.0, 600.0, 1330.0, 0.4, 2.5, 2.5, 2.25, 9e-7, 9e-7, 0.0, 125_000.0,
    )

    (
        relative_base_level,
        gridy_ref, temp_gridy_ref, density_gridy_ref, pressure_gridy_ref,
    ) = calculate_isostatic_relative_base_level(
        model_lith_thickness.total_column_thickness_meters,
        pressure_at_base_from_markers,
        reference_lith_thickness, reference_lith_thermal, density_model;
        iuse_linear_segments=iuse_linear_segments,
    )

    # Analytical isostatic check: with zero expansivity / compressibility and
    # equal continental-crust thicknesses (model 22+10 km, reference 22+10 km),
    # the analytical relative base level is identically zero.
    rho_w = 1000.0
    rho_c = density_model.upr_cont_crust.standard_density
    rho_m = density_model.asthenosphere.standard_density
    t_crust_model = (model_lith_thickness.thickness_upr_cont_crust_meters
                     + model_lith_thickness.thickness_lwr_cont_crust_meters)
    t_crust_ref = (reference_lith_thickness.thickness_upr_cont_crust_meters
                   + reference_lith_thickness.thickness_lwr_cont_crust_meters)
    relative_base_level_analytical = -(t_crust_ref - t_crust_model) *
                                      (rho_m - rho_c) / (rho_m - rho_w)

    expected_pressure_grid_base_pa    = 4.0471550000000005e9
    expected_pressure_markers_base_pa = 4.049276765625e9
    expected_pressure_ref_base_pa     = 4.370775500000001e9
    expected_relative_base_level_m    = -84.35073757759743

    @test isapprox(pressure_gridy[end],
                   expected_pressure_grid_base_pa; atol=1e-3)
    @test isapprox(pressure_at_base_from_markers,
                   expected_pressure_markers_base_pa; atol=1e-3)
    @test isapprox(pressure_gridy_ref[end],
                   expected_pressure_ref_base_pa; atol=1e-3)
    @test isapprox(relative_base_level,
                   expected_relative_base_level_m; atol=1e-6)

    # Sanity cross-check against the simple-crust analytical formula. The
    # formula only accounts for the crustal-thickness difference (zero here,
    # since both columns carry 32 km of crust). The numerical RBL also
    # absorbs the 10 km difference in total column thickness (model 130 km
    # vs reference 140 km of mantle), which the simple formula ignores —
    # so the two agree only to within ~100 m. Tolerance is loose by design;
    # this catches sign flips and order-of-magnitude regressions.
    @test isapprox(relative_base_level,
                   relative_base_level_analytical; atol=200.0)
end
