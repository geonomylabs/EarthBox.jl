using EarthBox
using Random
using Test

import EarthBox.Compaction: CompactionCorrection
import EarthBox.Compaction.MarkerCompaction: calculate_marker_swarm_indices
import EarthBox.Compaction.MarkerCompaction: calculate_x_sorted_swarm_indices_from_marker_x
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_layer_opt
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm_opt
import EarthBox.Markers.MarkerCoordinates.InitManager.MarkersRandomized: calc_coordinate_random

const _MCT_SEED = 20260501

const _MCT_X_DOMAIN          = 500_000.0
const _MCT_Y_DOMAIN          =  18_000.0
const _MCT_DX_MARKER         =     200.0   # scaled from demo's 50.0
const _MCT_DY_MARKER         =     200.0
const _MCT_STICKY_THICKNESS  =  10_000.0
const _MCT_WATER_THICKNESS   =   1_000.0
const _MCT_SEDIMENT_THICK    =   4_000.0
const _MCT_DX_TOPO           =     250.0
const _MCT_NEW_SED_THICK     =     200.0
const _MCT_POROSITY_INITIAL  =       0.5
const _MCT_DECAY_DEPTH       =   2_000.0
const _MCT_POROSITY_TRANSPORT = 0.5
const _MCT_DECAY_DEPTH_TRANSPORT = 2_000.0
const _MCT_NSMOOTH_TOP_BOTTOM = 2

const _MCT_MATID_STICKY      = Int16(1)
const _MCT_MATID_WATER       = Int16(2)
const _MCT_MATID_SED         = Int16(3)
const _MCT_MATID_MIDDLE_SED  = Int16(4)
const _MCT_MATID_CRUST       = Int16(5)

const _MCT_SED_BASIN_IDS = Int16[_MCT_MATID_SED, _MCT_MATID_MIDDLE_SED]
const _MCT_STICKY_IDS    = Int16[_MCT_MATID_STICKY, _MCT_MATID_WATER]

function _mct_define_marker_coords_random(nx::Int, ny::Int, dx::Float64, dy::Float64)
    nmarkers = nx * ny
    marker_x = zeros(nmarkers)
    marker_y = zeros(nmarkers)
    rns = rand(nmarkers)
    k = 1
    for j in 1:ny, i in 1:nx
        marker_x[k] = calc_coordinate_random(i - 1, dx, rns[k])
        marker_y[k] = calc_coordinate_random(j - 1, dy, rns[k])
        k += 1
    end
    return marker_x, marker_y
end

function _mct_make_markers(xmin_sed, xmax_sed)
    nx = floor(Int, _MCT_X_DOMAIN / _MCT_DX_MARKER) + 1
    ny = floor(Int, _MCT_Y_DOMAIN / _MCT_DY_MARKER) + 1
    nmarkers = nx * ny

    ymin_sed = _MCT_STICKY_THICKNESS + _MCT_WATER_THICKNESS
    ymax_sed = ymin_sed + _MCT_SEDIMENT_THICK
    ymin_mid = ymin_sed + 0.25 * _MCT_SEDIMENT_THICK
    ymax_mid = ymin_sed + 0.50 * _MCT_SEDIMENT_THICK

    marker_x, marker_y = _mct_define_marker_coords_random(
        nx, ny, _MCT_DX_MARKER, _MCT_DY_MARKER)

    matid = fill(_MCT_MATID_CRUST, nmarkers)
    porosity_initial = fill(Float32(_MCT_POROSITY_INITIAL), nmarkers)
    decay_depth = fill(Float32(_MCT_DECAY_DEPTH), nmarkers)
    max_burial_depth = zeros(Float32, nmarkers)

    for i in 1:nmarkers
        matid[i] = marker_y[i] <= _MCT_STICKY_THICKNESS ?
                   _MCT_MATID_STICKY : _MCT_MATID_CRUST
    end
    for i in 1:nmarkers
        if xmin_sed <= marker_x[i] <= xmax_sed &&
           _MCT_STICKY_THICKNESS < marker_y[i] < ymin_sed
            matid[i] = _MCT_MATID_WATER
        end
    end
    for i in 1:nmarkers
        if xmin_sed <= marker_x[i] <= xmax_sed &&
           ymin_sed <= marker_y[i] < ymax_sed
            matid[i] = _MCT_MATID_SED
        end
    end
    # Middle non-compacting layer (porosity_initial=0, decay_depth=1).
    for i in 1:nmarkers
        if ymin_mid < marker_y[i] < ymax_mid && matid[i] == _MCT_MATID_SED
            matid[i] = _MCT_MATID_MIDDLE_SED
            porosity_initial[i] = Float32(0.0)
            decay_depth[i] = Float32(1.0)
        end
    end

    return marker_x, marker_y, matid, porosity_initial, decay_depth, max_burial_depth
end

function _mct_setup(xmin_sed, xmax_sed)
    marker_x, marker_y, matid, porosity_initial, decay_depth, max_burial_depth =
        _mct_make_markers(xmin_sed, xmax_sed)

    topo_gridx = collect(0.0:_MCT_DX_TOPO:_MCT_X_DOMAIN - _MCT_DX_TOPO)
    topo_gridy = fill(_MCT_STICKY_THICKNESS, length(topo_gridx))
    for i in eachindex(topo_gridx)
        if xmin_sed <= topo_gridx[i] <= xmax_sed
            topo_gridy[i] = _MCT_STICKY_THICKNESS + _MCT_WATER_THICKNESS
        end
    end

    new_sed_thickness_gridx = fill(_MCT_NEW_SED_THICK, length(topo_gridx))
    for i in eachindex(topo_gridx)
        if topo_gridx[i] < xmin_sed || topo_gridx[i] >= xmax_sed
            new_sed_thickness_gridx[i] = 0.0
        end
    end

    indices_sed_basin = calculate_marker_swarm_indices(marker_x, matid, _MCT_SED_BASIN_IDS)
    indices_sticky    = calculate_marker_swarm_indices(marker_x, matid, _MCT_STICKY_IDS)

    x_sorted_indices_sticky = calculate_x_sorted_swarm_indices_from_marker_x(
        marker_x, indices_sticky)
    x_sorted_indices_sed_basin = calculate_x_sorted_swarm_indices_from_marker_x(
        marker_x, indices_sed_basin)

    search_radius = _MCT_DX_TOPO / 2.0

    top_sed, bottom_sed = calculate_top_and_bottom_of_layer_opt(
        _MCT_SED_BASIN_IDS, matid, marker_x, marker_y, topo_gridx, search_radius)
    sed_thickness_gridx = bottom_sed .- top_sed

    top_sticky, bottom_sticky = calculate_top_and_bottom_of_swarm_opt(
        x_sorted_indices_sticky, marker_x, marker_y, topo_gridx, search_radius)
    sticky_thickness_gridx = bottom_sticky .- top_sticky

    return (
        marker_x, marker_y, matid,
        porosity_initial, decay_depth, max_burial_depth,
        topo_gridx, topo_gridy, new_sed_thickness_gridx,
        indices_sed_basin, indices_sticky, x_sorted_indices_sed_basin,
        sed_thickness_gridx, sticky_thickness_gridx,
    )
end

@testset "MarkerCompaction (3-step correction)" begin
    Random.seed!(_MCT_SEED)

    xmin_sed = _MCT_X_DOMAIN * 0.25
    xmax_sed = _MCT_X_DOMAIN * 0.75

    (
        marker_x, marker_y, matid,
        porosity_initial, decay_depth, max_burial_depth,
        topo_gridx, topo_gridy, new_sed_thickness_gridx,
        indices_sed_basin, indices_sticky, x_sorted_indices_sed_basin,
        sed_thickness_gridx, sticky_thickness_gridx,
    ) = _mct_setup(xmin_sed, xmax_sed)

    nmarkers = length(marker_x)
    ntopo = length(topo_gridx)
    compaction_array_buf = zeros(Float64, ntopo, 20, 9)
    markers_topo_xindex_buf = Vector{Int64}(undef, nmarkers)
    markers_compaction_yindex_buf = Vector{Int64}(undef, nmarkers)
    markers_unit_distance_from_cell_top_buf = Vector{Float64}(undef, nmarkers)
    total_marker_compaction_displacement_buf = Vector{Float64}(undef, nmarkers)
    sticky_displacement_factors_buf = Vector{Float64}(undef, nmarkers)
    sticky_marker_displacement_buf = Vector{Float64}(undef, nmarkers)

    topo_gridy_initial = copy(topo_gridy)
    sediment_thickness_initial = copy(sed_thickness_gridx)
    sticky_thickness_initial = copy(sticky_thickness_gridx)
    topo_gridy_transport = topo_gridy .- new_sed_thickness_gridx

    porosity_initial_f32 = Float32.(porosity_initial)
    decay_depth_f32 = Float32.(decay_depth)
    max_burial_depth_f32 = Float32.(max_burial_depth)

    search_radius = _MCT_DX_TOPO / 2.0

    local topo_gridy_corrected, total_sed_thickness_corrected, new_thickness_decompacted

    for istep in 1:3
        topo_gridy_corrected, total_sed_thickness_corrected, new_thickness_decompacted =
            CompactionCorrection.decompact_new_sediment_and_compact_markers(
                compaction_array_buf,
                markers_topo_xindex_buf,
                markers_compaction_yindex_buf,
                markers_unit_distance_from_cell_top_buf,
                total_marker_compaction_displacement_buf,
                sticky_displacement_factors_buf,
                sticky_marker_displacement_buf,
                _MCT_POROSITY_TRANSPORT, _MCT_DECAY_DEPTH_TRANSPORT,
                topo_gridx, topo_gridy_initial, topo_gridy_transport,
                marker_x, marker_y,
                porosity_initial_f32, decay_depth_f32, max_burial_depth_f32,
                sediment_thickness_initial, sticky_thickness_initial,
                indices_sed_basin, indices_sticky,
                x_sorted_indices_sed_basin,
                search_radius,
            )

        # Update state for next step (mirrors demo).
        topo_gridy_initial = copy(topo_gridy_corrected)
        topo_gridy_transport = topo_gridy_initial .- new_sed_thickness_gridx
        sediment_thickness_initial = copy(total_sed_thickness_corrected)
        # Recompute sticky thickness from current marker positions.
        top_sticky, bottom_sticky = calculate_top_and_bottom_of_layer_opt(
            _MCT_STICKY_IDS, matid, marker_x, marker_y, topo_gridx, search_radius)
        sticky_thickness_initial = bottom_sticky .- top_sticky
    end

    @test length(topo_gridy_corrected) == ntopo
    @test length(total_sed_thickness_corrected) == ntopo
    @test length(new_thickness_decompacted) == ntopo

    # Sample indices spanning sticky-only → basin → sticky-only.
    # x_domain = 500 km, dx_topo = 250 m → ntopo = 2000. Basin in [125 km, 375 km]
    # ⇒ topo indices ~501..1500.
    sample_indices = [1, 501, 1001, 1500, 2000]

    expected_topo_gridy_corrected = [
        10000.0,
        10338.58449611313,
        10286.669809950608,
        10251.158914305868,
        10000.0,
    ]
    expected_total_sed_thickness_corrected = [
        0.0,
        3604.5703819778637,
        4526.877108222951,
        2952.295672962172,
        0.0,
    ]

    for (k, idx) in enumerate(sample_indices)
        @test isapprox(topo_gridy_corrected[idx],
                       expected_topo_gridy_corrected[k]; atol=1e-6)
        @test isapprox(total_sed_thickness_corrected[idx],
                       expected_total_sed_thickness_corrected[k]; atol=1e-6)
    end

    expected_min_topo_corrected         = 10000.0
    expected_max_topo_corrected         = 11000.0
    expected_min_thickness_decompacted  = 0.0
    expected_max_thickness_decompacted  = 365.76919209877974
    expected_marker_y_sed_checksum      = 3.268929086338946e8

    @test isapprox(minimum(topo_gridy_corrected),
                   expected_min_topo_corrected; atol=1e-6)
    @test isapprox(maximum(topo_gridy_corrected),
                   expected_max_topo_corrected; atol=1e-6)
    @test isapprox(minimum(new_thickness_decompacted),
                   expected_min_thickness_decompacted; atol=1e-6)
    @test isapprox(maximum(new_thickness_decompacted),
                   expected_max_thickness_decompacted; atol=1e-6)

    # Marker advection checksum: catches regressions in marker displacement
    # that wouldn't show up in topography-only assertions.
    marker_y_sed_checksum = sum(marker_y[indices_sed_basin])
    @test isapprox(marker_y_sed_checksum,
                   expected_marker_y_sed_checksum; atol=1e-3)
end
