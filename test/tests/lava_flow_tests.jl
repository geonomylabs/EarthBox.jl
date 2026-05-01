using EarthBox
using Random
using Test

import EarthBox.ModelDataContainer: ModelData
import EarthBox.SurfaceProcesses.LavaFlowManager: run_lava_flow_model

const _LFT_SEED = 20260429

const _LFT_XSIZE   = 500_000.0
const _LFT_YSIZE   =  30_000.0
const _LFT_XNUM    = 21
const _LFT_YNUM    = 11
const _LFT_TOPONUM = 1001
const _LFT_DX_TOPO = _LFT_XSIZE / Float64(_LFT_TOPONUM - 1)

const _LFT_X_BASIN_MIN = 100_000.0
const _LFT_X_BASIN_MAX = 400_000.0

const _LFT_Y_SEALEVEL          = 5_000.0
const _LFT_STICKY_THICKNESS    = 5_000.0
const _LFT_SEDIMENT_THICKNESS  = 4_000.0
const _LFT_BASEMENT_DEPTH      = _LFT_Y_SEALEVEL + _LFT_STICKY_THICKNESS + _LFT_SEDIMENT_THICKNESS

const _LFT_ERUPTION_X_MIN              = 220_000.0
const _LFT_ERUPTION_DOMAIN_WIDTH       =  60_000.0
const _LFT_TOTAL_EXTRUSION_VOLUME      = 2.0 * 48_000_000.0
const _LFT_RESIDUAL_THICKNESS          =  30.0
const _LFT_CHARACTERISTIC_FLOW_LENGTH  =  60_000.0
const _LFT_DECIMATION_FACTOR           = 4

const _LFT_SEDIMENT_POROSITY_INITIAL = 0.7
const _LFT_SEDIMENT_DECAY_DEPTH      = 2_500.0
const _LFT_SEDIMENT_PRECOMPACTED     = true

const _LFT_MARKER_NX = _LFT_TOPONUM
const _LFT_MARKER_NY = ceil(Int, _LFT_YSIZE / 100.0) + 1

function _lft_populate_material_types!(model::ModelData)
    types = model.materials.dicts.matid_types
    types["StickyAir"]                              = Int16[]
    types["StickyWater"]                            = Int16[Int16(2)]
    types["Sediment"]                               = Int16[Int16(3)]
    types["SolidifiedBasalt"]                       = Int16[]
    types["Salt"]                                   = Int16[]
    types["ExtractedGabbroicMagma"]                 = Int16[Int16(91)]
    types["SolidifiedGabbroPartiallyMolten"]        = Int16[Int16(92)]
    types["ExtractedLayeredGabbroicMagma"]          = Int16[Int16(93)]
    types["SolidifiedLayeredGabbroPartiallyMolten"] = Int16[Int16(94)]
    return nothing
end

function _lft_configure_extrusion_parameters!(model::ModelData)
    extrusion = model.melting.parameters.extrusion
    extrusion.iuse_extrusion.value                            = 1
    extrusion.iuse_random_eruption_location.value             = 1
    extrusion.iuse_normal_eruption_location.value             = 0
    extrusion.decimation_factor.value                         = _LFT_DECIMATION_FACTOR
    extrusion.characteristic_flow_length_subaerial.value      = _LFT_CHARACTERISTIC_FLOW_LENGTH
    extrusion.characteristic_flow_length_submarine.value      = _LFT_CHARACTERISTIC_FLOW_LENGTH
    extrusion.residual_lava_thickness_subaerial.value         = _LFT_RESIDUAL_THICKNESS
    extrusion.residual_lava_thickness_submarine.value         = _LFT_RESIDUAL_THICKNESS
    extrusion.porosity_initial_lava_flow.value                = 0.4
    extrusion.decay_depth_lava_flow.value                     = 2_000.0
    extrusion.initial_magma_flush_steps.value                 = 0

    extrusion.width_eruption_domain_fixed.value     = _LFT_ERUPTION_DOMAIN_WIDTH
    extrusion.width_eruption_domain_fixed_max.value = _LFT_ERUPTION_DOMAIN_WIDTH
    extrusion.characteristic_magmatic_crust_height.value     = 1.0
    extrusion.characteristic_magmatic_crust_height_min.value = 0.0
    extrusion.characteristic_magmatic_crust_height_max.value = 2.0

    model.topography.parameters.sealevel.y_sealevel.value = _LFT_Y_SEALEVEL
    return nothing
end

function _lft_configure_compaction_parameters!(model::ModelData)
    model.topography.parameters.downhill_diffusion.iuse_compaction_correction.value = 1
    return nothing
end

function _lft_initialize_topography!(model::ModelData)
    gridt = model.topography.arrays.gridt.array
    toponum = size(gridt, 2)
    @assert toponum == _LFT_TOPONUM
    for i in 1:toponum
        x = (i - 1) * _LFT_DX_TOPO
        gridt[1, i] = x
        gridt[2, i] = _LFT_Y_SEALEVEL + _LFT_STICKY_THICKNESS
        gridt[7, i] = 0.0
    end
    return nothing
end

function _lft_populate_markers!(model::ModelData)
    matid_types = model.materials.dicts.matid_types
    matid_sticky_water = matid_types["StickyWater"][1]
    matid_sediment     = matid_types["Sediment"][1]

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_porosity_initial = model.markers.arrays.material.marker_porosity_initial.array
    marker_decay_depth      = model.markers.arrays.material.marker_decay_depth.array
    marker_max_burial_depth = model.markers.arrays.material.marker_max_burial_depth.array

    nmarkers = length(marker_x)
    nx_m = _LFT_MARKER_NX
    ny_m = _LFT_MARKER_NY
    nm = nx_m * ny_m
    @assert nm <= nmarkers
    dx_m = _LFT_XSIZE / Float64(nx_m)
    dy_m = _LFT_YSIZE / Float64(ny_m)

    for i in 1:nmarkers
        marker_x[i] = 0.0
        marker_y[i] = 0.0
        marker_matid[i] = Int16(0)
        marker_porosity_initial[i] = Float32(0.0)
        marker_decay_depth[i]      = Float32(2_000.0)
        marker_max_burial_depth[i] = Float32(0.0)
    end

    k = 0
    for j in 1:ny_m, i in 1:nx_m
        k += 1
        if k > nm
            break
        end
        x = (i - 0.5) * dx_m
        y = (j - 0.5) * dy_m
        marker_x[k] = x
        marker_y[k] = y
        if y < _LFT_Y_SEALEVEL + _LFT_STICKY_THICKNESS
            marker_matid[k] = matid_sticky_water
            marker_porosity_initial[k] = Float32(0.0)
            marker_decay_depth[k] = Float32(2_000.0)
            marker_max_burial_depth[k] = Float32(0.0)
        elseif x >= _LFT_X_BASIN_MIN && x <= _LFT_X_BASIN_MAX && y <= _LFT_BASEMENT_DEPTH
            marker_matid[k] = matid_sediment
            marker_porosity_initial[k] = Float32(_LFT_SEDIMENT_POROSITY_INITIAL)
            marker_decay_depth[k] = Float32(_LFT_SEDIMENT_DECAY_DEPTH)
            marker_max_burial_depth[k] = _LFT_SEDIMENT_PRECOMPACTED ?
                Float32(y - (_LFT_Y_SEALEVEL + _LFT_STICKY_THICKNESS)) : Float32(0.0)
        else
            marker_matid[k] = Int16(0)
            marker_porosity_initial[k] = Float32(0.0)
            marker_decay_depth[k] = Float32(2_000.0)
            marker_max_burial_depth[k] = Float32(0.0)
        end
    end
    return nothing
end

function _lft_configure_drainage_basin!(model::ModelData)
    model.melting.parameters.extraction.ndrainage_basin.value = 1
    arrays = model.melting.arrays.extraction
    arrays.xstart_drainage.array[1]                  = 0.0
    arrays.xend_drainage.array[1]                    = _LFT_XSIZE
    arrays.extrusion_volumes.array[1]                = _LFT_TOTAL_EXTRUSION_VOLUME
    arrays.xmid_molten_zones.array[1]                = _LFT_ERUPTION_X_MIN + _LFT_ERUPTION_DOMAIN_WIDTH / 2.0
    arrays.width_molten_zones.array[1]               = _LFT_ERUPTION_DOMAIN_WIDTH
    arrays.avg_shallow_partial_melt_xcoors.array[1]  = _LFT_ERUPTION_X_MIN + _LFT_ERUPTION_DOMAIN_WIDTH / 2.0
    arrays.avg_shallow_partial_melt_ycoors.array[1]  = _LFT_Y_SEALEVEL + _LFT_STICKY_THICKNESS
    return nothing
end

function _lft_build_model()::ModelData
    init_params = Dict{String, Vector{Any}}(
        "xnum"             => Any[_LFT_XNUM,  "None"],
        "ynum"             => Any[_LFT_YNUM,  "None"],
        "xsize"            => Any[_LFT_XSIZE, "m"],
        "ysize"            => Any[_LFT_YSIZE, "m"],
        "nmarkers_cell_x"  => Any[50.0,       "None"],
        "nmarkers_cell_y"  => Any[35.0,       "None"],
    )
    model = ModelData(nothing, init_params)
    _lft_populate_material_types!(model)
    _lft_configure_extrusion_parameters!(model)
    _lft_configure_compaction_parameters!(model)
    _lft_initialize_topography!(model)
    _lft_populate_markers!(model)
    _lft_configure_drainage_basin!(model)
    return model
end

function _lft_build_sediment_thickness_initial()::Vector{Float64}
    sediment_thickness = zeros(Float64, _LFT_TOPONUM)
    for i in 1:_LFT_TOPONUM
        x = (i - 1) * _LFT_DX_TOPO
        if x >= _LFT_X_BASIN_MIN && x <= _LFT_X_BASIN_MAX
            sediment_thickness[i] = _LFT_SEDIMENT_THICKNESS
        end
    end
    return sediment_thickness
end

@testset "LavaFlowCompactionThickSediment" begin
    Random.seed!(_LFT_SEED)
    model = _lft_build_model()
    sediment_thickness_initial = _lft_build_sediment_thickness_initial()

    gridt = model.topography.arrays.gridt.array
    pre_topo_y = copy(gridt[2, :])

    run_lava_flow_model(model, sediment_thickness_initial)

    extrusion_thickness = copy(gridt[7, :])
    post_topo_y = copy(gridt[2, :])
    topo_delta = post_topo_y .- pre_topo_y

    @test length(extrusion_thickness) == _LFT_TOPONUM

    sample_indices = [1, 401, 451, 501, 551, 601, 1001]

    expected_extrusion_thickness = [
        0.0,
        744.5189355060892,
        1168.7585749596215,
        1161.6746656890068,
        1084.5360997270568,
        804.5738524235875,
        0.0,
    ]
    expected_max_extrusion_thickness  = 1175.2321311439564
    expected_total_lava_volume        = 1.4040035336284637e8
    expected_nonzero_extrusion_count  = 335
    expected_max_topo_subsidence      = 0.0
    expected_min_topo_subsidence      = -595.606651107064

    for (k, idx) in enumerate(sample_indices)
        @test isapprox(extrusion_thickness[idx],
                       expected_extrusion_thickness[k]; atol=1e-6)
    end
    @test isapprox(maximum(extrusion_thickness),
                   expected_max_extrusion_thickness; atol=1e-6)
    @test isapprox(sum(extrusion_thickness) * _LFT_DX_TOPO,
                   expected_total_lava_volume; atol=1e-3)
    @test count(>(0.0), extrusion_thickness) == expected_nonzero_extrusion_count
    @test isapprox(maximum(topo_delta), expected_max_topo_subsidence; atol=1e-6)
    @test isapprox(minimum(topo_delta), expected_min_topo_subsidence; atol=1e-6)
end
