using EarthBox
using Test

import EarthBox.Gravity: gravity_anomaly_loop_left_edge, calculate_free_air_gravity
import EarthBox.Gravity.GravityCell: calc_cell_gravity_for_xarray, calculate_gravity_anomaly_prism

# Inputs from Fowler 2001, Figure 5.6a (Bott 1982): mountain over an
# isostatically-compensated crustal root.
const _GT_XSIZE       = 700_000.0
const _GT_YSIZE       = 210_000.0
const _GT_DX          = 1_000.0
const _GT_DY          = 1_000.0
const _GT_THICK_AIR   =  10_000.0
const _GT_THICK_CRUST =  30_000.0
const _GT_RHO_CRUST   = 2850.0
const _GT_RHO_MANTLE  = 3300.0
const _GT_RHO_AIR     =    0.0
const _GT_RHO_WATER   = 1000.0
const _GT_HEIGHT_TOPO =  3_000.0
const _GT_XO_TOPO     = 100_000.0
const _GT_WIDTH_TOPO  = 500_000.0
const _GT_GRAVITY_CONSTANT = 6.6732e-11

function _gt_build_grids()
    xnum = floor(Int, _GT_XSIZE / _GT_DX) + 1
    ynum = floor(Int, _GT_YSIZE / _GT_DY) + 1

    bgridx = collect(range(0.0, _GT_XSIZE, length=xnum))
    bgridy = collect(range(0.0, _GT_YSIZE, length=ynum))
    xsteps = diff(bgridx)
    ysteps = diff(bgridy)

    rho_grid = fill(_GT_RHO_MANTLE, ynum, xnum)
    topo_gridx = collect(range(0.0, _GT_XSIZE, length=xnum))
    topo_gridy = fill(_GT_THICK_AIR, xnum)

    # Mountain topography (negative-up convention: y is depth-positive).
    for i in 1:xnum
        if _GT_XO_TOPO <= topo_gridx[i] <= _GT_XO_TOPO + _GT_WIDTH_TOPO
            topo_gridy[i] -= _GT_HEIGHT_TOPO
        end
    end

    # Air layer (everything above sea level).
    for j in 1:xnum, i in 1:ynum
        if bgridy[i] < _GT_THICK_AIR
            rho_grid[i, j] = _GT_RHO_AIR
        end
    end

    # Initial flat crust (sea level → moho).
    y_moho = _GT_THICK_AIR + _GT_THICK_CRUST
    for j in 1:xnum, i in 1:ynum
        if _GT_THICK_AIR <= bgridy[i] <= y_moho
            rho_grid[i, j] = _GT_RHO_CRUST
        end
    end

    # Mountain + isostatic root.
    root_thickness = _GT_HEIGHT_TOPO * _GT_RHO_CRUST /
                     (_GT_RHO_MANTLE - _GT_RHO_CRUST)
    y_topo_top = _GT_THICK_AIR - _GT_HEIGHT_TOPO
    y_root_bot = _GT_THICK_AIR + _GT_THICK_CRUST + root_thickness
    for j in 1:xnum, i in 1:ynum
        x_grid = bgridx[j]
        y_grid = bgridy[i]
        if _GT_XO_TOPO <= x_grid <= _GT_XO_TOPO + _GT_WIDTH_TOPO
            if y_topo_top <= y_grid <= y_root_bot
                rho_grid[i, j] = _GT_RHO_CRUST
            end
        end
    end

    return bgridx, bgridy, xsteps, ysteps, rho_grid, topo_gridx, topo_gridy
end

@testset "GridGravity (mountain + root)" begin
    bgridx, bgridy, xsteps, ysteps, rho_grid, topo_gridx, topo_gridy = _gt_build_grids()

    bouguer_mgal = gravity_anomaly_loop_left_edge(
        _GT_THICK_AIR, bgridx, bgridy, xsteps, ysteps, rho_grid,
        topo_gridx, _GT_GRAVITY_CONSTANT,
    )
    free_air_mgal = calculate_free_air_gravity(
        _GT_THICK_AIR, _GT_GRAVITY_CONSTANT, _GT_RHO_CRUST,
        topo_gridx, topo_gridy, bouguer_mgal,
    )

    @test length(bouguer_mgal) == length(topo_gridx)
    @test length(free_air_mgal) == length(topo_gridx)

    # Sample indices spanning off-mountain → mountain interior → off-mountain.
    # Topography lives between x=100 km and x=600 km; with dx=1 km, indices
    # 101..601 are over the mountain.
    sample_indices = [1, 101, 201, 351, 501, 601, 701]

    expected_bouguer_mgal = [
        -35.92560434193713,
        -171.60095487337216,
        -304.01151186331737,
        -322.3713614309996,
        -304.01151186331094,
        -171.60095487337415,
        -35.92560434192124,
    ]
    expected_free_air_mgal = [
        -35.92560434193713,
         186.89158636712332,
          54.48102937717812,
          36.12117980949586,
          54.48102937718454,
         186.89158636712133,
         -35.92560434192124,
    ]
    expected_min_bouguer  = -322.3713614309996
    expected_max_bouguer  =  -35.92560434192124
    expected_min_free_air = -168.7108619452204
    expected_max_free_air =  186.89158636712332
    expected_sum_bouguer  = -163470.00613723617

    for (k, idx) in enumerate(sample_indices)
        @test isapprox(bouguer_mgal[idx],  expected_bouguer_mgal[k];  atol=1e-6)
        @test isapprox(free_air_mgal[idx], expected_free_air_mgal[k]; atol=1e-6)
    end
    @test isapprox(minimum(bouguer_mgal),  expected_min_bouguer;  atol=1e-6)
    @test isapprox(maximum(bouguer_mgal),  expected_max_bouguer;  atol=1e-6)
    @test isapprox(minimum(free_air_mgal), expected_min_free_air; atol=1e-6)
    @test isapprox(maximum(free_air_mgal), expected_max_free_air; atol=1e-6)
    @test isapprox(sum(bouguer_mgal),      expected_sum_bouguer;  atol=1e-3)
end

# Cross-check the finite-prism kernel (Turcotte & Schubert 2014) against the
# infinite-beam analytical solution (Telford et al. 1990) for a thin, very
# long y-extruded prism. Because the box is 200 m wide in y but only 1 m
# wide in x, the y-extrusion is large enough that the finite prism should
# match the infinite-beam closely. This is a correctness test, not a
# regression pin.
@testset "CellGravityPrism (Telford vs Turcotte)" begin
    G = 6.6732e-11
    delta_density = 1.0

    xmin_observer = -3.0
    xmax_observer =  3.0
    xspacing_observer = 0.1
    nobs = floor(Int, (xmax_observer - xmin_observer) / xspacing_observer) + 1
    xcoors_observer = collect(xmin_observer:xspacing_observer:xmax_observer)
    @assert length(xcoors_observer) == nobs

    xwidth_box = 1.0
    x1_box = -xwidth_box / 2.0
    x2_box =  xwidth_box / 2.0
    y1_box = -100.0
    y2_box =  100.0
    z1_box = 1.0 / 3.0
    z2_box = 4.0 / 3.0

    # Telford infinite-beam (extruded in y) solution.
    z_datum = 0.0
    cell_width = x2_box - x1_box
    cell_height = z2_box - z1_box
    grav_anomaly_x = calc_cell_gravity_for_xarray(
        G, delta_density,
        x1_box, z1_box, z_datum,
        cell_height, cell_width,
        xcoors_observer,
    )
    nondim_telford = grav_anomaly_x ./ (2.0 * G * delta_density)

    # Turcotte & Schubert finite-prism solution at every observer.
    nondim_turcotte = zeros(nobs)
    ycoor_observer = 0.0
    zcoor_observer = 0.0
    for i in 1:nobs
        x = xcoors_observer[i]
        dx1 = x1_box - x;  dx2 = x2_box - x
        dy1 = y1_box - ycoor_observer;  dy2 = y2_box - ycoor_observer
        dz1 = z1_box - zcoor_observer;  dz2 = z2_box - zcoor_observer
        delta_g = calculate_gravity_anomaly_prism(
            G, delta_density, dx1, dx2, dy1, dy2, dz1, dz2,
        )
        nondim_turcotte[i] = delta_g / (2.0 * G * delta_density)
    end

    # Cross-check: the two solutions should agree to better than 1% at every
    # observer because the prism is very long in y relative to its x extent.
    @test length(nondim_telford) == length(nondim_turcotte) == nobs
    for i in 1:nobs
        @test isapprox(nondim_turcotte[i], nondim_telford[i]; rtol=1e-2, atol=1e-6)
    end
end
