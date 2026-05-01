"""
Integration test: run_lava_flow_model on a thick pre-existing sediment basin
with compaction correction enabled.

This test exercises the lava-flow path that none of the other lava-flow tests
reach: extrude_magma → apply_compaction_model! → marker compaction with
non-trivial sediment thickness. It builds a real ModelData (the only one of
the lava-flow tests to do so), populates a thick sediment basin via markers
and topography, runs run_lava_flow_model for one synthetic eruption event,
and emits before/after plots for visual validation.

There are no bit-equivalence assertions — visual inspection of the output
PNGs is the verification.

Outputs:
- compaction_thick_sediment_topography.png — topography before vs after.
- compaction_thick_sediment_lava.png — total lava thickness vs x.
- compaction_thick_sediment_markers.png — marker positions before vs after.
"""
module LavaFlowTestCompactionThickSediment

using CairoMakie
import EarthBox.ModelDataContainer: ModelData
import EarthBox.SurfaceProcesses.LavaFlowManager: run_lava_flow_model

const SEED = 20260429

# ── Geometry & resolution ───────────────────────────────────────────────────
const XSIZE   = 500_000.0      # m, domain x-extent
const YSIZE   =  30_000.0      # m, domain y-extent
const XNUM    = 21             # basic grid x-nodes
const YNUM    = 11             # basic grid y-nodes

# Topography grid (lava-flow grid). The default toponum from the parameter
# registry is 1001 — derive DX_TOPO from that.
const TOPONUM = 1001
const DX_TOPO = XSIZE / Float64(TOPONUM - 1)

# Sediment basin extent in x.
const X_BASIN_MIN = 100_000.0
const X_BASIN_MAX = 400_000.0

# Vertical layout (y is depth; positive y is below sea level).
const Y_SEALEVEL          = 5_000.0   # m below model y=0
const STICKY_THICKNESS    = 5_000.0   # sticky water/air column above seabed
const SEDIMENT_THICKNESS  =  4_000.0   # the basin fill we're testing on
const BASEMENT_DEPTH      = Y_SEALEVEL + STICKY_THICKNESS + SEDIMENT_THICKNESS  # 24_000 m

# ── Eruption config ─────────────────────────────────────────────────────────
const ERUPTION_X_MIN              = 220_000.0
const ERUPTION_DOMAIN_WIDTH       =  60_000.0
const TOTAL_EXTRUSION_VOLUME      = 2.0*48_000_000.0   # m^2 (per-unit-z volume)
const RESIDUAL_THICKNESS          =  30.0
const CHARACTERISTIC_FLOW_LENGTH  =  60_000.0
const DECIMATION_FACTOR           = 4

# ── Sediment compaction (Athy's law: φ(z) = φ₀ · exp(−z / λ)) ───────────────
# Bump φ₀ up and λ down for a squishier basin. Set SEDIMENT_PRECOMPACTED=false
# to start the basin "fresh" (no prior burial credit) so even small lava loads
# trigger compaction; true keeps the original behavior where each marker's
# current depth counts as already-compacted.
const SEDIMENT_POROSITY_INITIAL   = 0.7      # φ₀ (fraction). Try 0.7–0.8 for squishy.
const SEDIMENT_DECAY_DEPTH        = 2_500.0  # λ (m).         Try 500–1_000 for squishy.
const SEDIMENT_PRECOMPACTED       = true     # false = max_burial_depth starts at 0.

# ── Marker layout ───────────────────────────────────────────────────────────
# The compaction grid has 20 cells per topography column (hardcoded in
# MarkerCompaction). For a 4 km basin those cells are 200 m tall. Empty cells
# fall back to the lava-flow decompaction porosity (extrusion.porosity_initial_lava_flow,
# 0.4), so we pack markers densely enough that every (topo column × compaction
# cell) contains at least one marker — this lets the per-marker SEDIMENT_*
# values actually drive the compaction kernel.
#   MARKER_NX: one marker per topo column (DX_TOPO spacing in x).
#   MARKER_NY: ~2 markers per 200 m compaction cell (≈100 m spacing in y).
const MARKER_NX = TOPONUM                            # 1001
const MARKER_NY = ceil(Int, YSIZE / 100.0) + 1        # 301

function run_test()
    println("=== LavaFlowTestCompactionThickSediment ===")
    model = build_model()
    sediment_thickness_initial = build_sediment_thickness_initial(model)

    println(">> Capturing pre-extrusion state...")
    pre = snapshot_state(model)

    println(">> Calling run_lava_flow_model...")
    run_lava_flow_model(model, sediment_thickness_initial)

    println(">> Capturing post-extrusion state...")
    post = snapshot_state(model)

    print_diagnostics(pre, post)

    println(">> Plotting...")
    plot_topography(pre, post)
    plot_lava_thickness(post)
    plot_markers(model, pre, post)

    println("Done. PNG outputs written to working directory.")
    return nothing
end

# ── ModelData construction ──────────────────────────────────────────────────

function build_model()::ModelData
    # nmarkers_cell_* sized so total marker slots ≥ MARKER_NX × MARKER_NY
    # (which is what populate_markers! actually uses). Total = nmcx × nmcy ×
    # (xnum-1) × (ynum-1) = 50 × 35 × 20 × 10 = 350_000 slots, comfortably
    # above the ~301k we lay out below.
    init_params = Dict{String, Vector{Any}}(
        "xnum"             => Any[XNUM,  "None"],
        "ynum"             => Any[YNUM,  "None"],
        "xsize"            => Any[XSIZE, "m"],
        "ysize"            => Any[YSIZE, "m"],
        "nmarkers_cell_x"  => Any[50.0,  "None"],
        "nmarkers_cell_y"  => Any[35.0,  "None"],
    )
    model = ModelData(nothing, init_params)

    populate_material_types!(model)
    configure_extrusion_parameters!(model)
    configure_compaction_parameters!(model)
    initialize_topography!(model)
    populate_markers!(model)
    configure_drainage_basin!(model)

    return model
end

# matid_types is empty by default — the production code populates it from the
# material registry. Populate the keys the lava-flow + compaction code paths
# look up (see MarkerCompaction.get_sticky_material_ids,
# get_sedimentary_basin_material_ids, MaterialGroupIDs.get_molten_gabbro_ids).
# Keys that production never indexes from this test path get empty vectors.
function populate_material_types!(model::ModelData)
    types = model.materials.dicts.matid_types
    types["StickyAir"]                              = Int16[]
    types["StickyWater"]                            = Int16[Int16(2)]
    types["Sediment"]                               = Int16[Int16(3)]
    types["SolidifiedBasalt"]                       = Int16[]
    types["Salt"]                                   = Int16[]
    # Gabbro keys are looked up by get_extrusion_location_parameters via
    # get_molten_gabbro_ids and calculate_dimensions_of_gabbro_glacier_domain
    # even when width_eruption_domain_fixed is set. Populate with placeholder
    # ids that no marker carries, so the find_shallowest path is a no-op.
    types["ExtractedGabbroicMagma"]                 = Int16[Int16(91)]
    types["SolidifiedGabbroPartiallyMolten"]        = Int16[Int16(92)]
    types["ExtractedLayeredGabbroicMagma"]          = Int16[Int16(93)]
    types["SolidifiedLayeredGabbroPartiallyMolten"] = Int16[Int16(94)]
    return nothing
end

function configure_extrusion_parameters!(model::ModelData)
    extrusion = model.melting.parameters.extrusion
    extrusion.iuse_extrusion.value                            = 1
    extrusion.iuse_random_eruption_location.value             = 1
    extrusion.iuse_normal_eruption_location.value             = 0
    extrusion.decimation_factor.value                         = DECIMATION_FACTOR
    extrusion.characteristic_flow_length_subaerial.value      = CHARACTERISTIC_FLOW_LENGTH
    extrusion.characteristic_flow_length_submarine.value      = CHARACTERISTIC_FLOW_LENGTH
    extrusion.residual_lava_thickness_subaerial.value         = RESIDUAL_THICKNESS
    extrusion.residual_lava_thickness_submarine.value         = RESIDUAL_THICKNESS
    extrusion.porosity_initial_lava_flow.value                = 0.4
    extrusion.decay_depth_lava_flow.value                     = 2_000.0
    extrusion.initial_magma_flush_steps.value                 = 0

    # Drive width selection through the fixed-width branch so we don't depend
    # on gabbro-glacier markers being present.
    extrusion.width_eruption_domain_fixed.value     = ERUPTION_DOMAIN_WIDTH
    extrusion.width_eruption_domain_fixed_max.value = ERUPTION_DOMAIN_WIDTH
    extrusion.characteristic_magmatic_crust_height.value     = 1.0
    extrusion.characteristic_magmatic_crust_height_min.value = 0.0
    extrusion.characteristic_magmatic_crust_height_max.value = 2.0

    model.topography.parameters.sealevel.y_sealevel.value = Y_SEALEVEL
    return nothing
end

function configure_compaction_parameters!(model::ModelData)
    model.topography.parameters.downhill_diffusion.iuse_compaction_correction.value = 1
    return nothing
end

function initialize_topography!(model::ModelData)
    gridt = model.topography.arrays.gridt.array
    toponum = size(gridt, 2)
    @assert toponum == TOPONUM "Topography grid size mismatch: $toponum vs $TOPONUM"

    for i in 1:toponum
        x = (i - 1) * DX_TOPO
        gridt[1, i] = x
        # Top of sediment = sealevel + sticky water column thickness; everywhere flat.
        gridt[2, i] = Y_SEALEVEL + STICKY_THICKNESS
        gridt[7, i] = 0.0
    end
    return nothing
end

function populate_markers!(model::ModelData)
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
    println(">> Marker array length: $nmarkers")

    # Lay out markers on a structured grid spanning the model x/y domain. The
    # sediment basin lives in [X_BASIN_MIN, X_BASIN_MAX] × [Y_SEALEVEL+STICKY,
    # BASEMENT_DEPTH]. Everything above the seabed is sticky water; basement
    # below the sediment is unmarked (left as the constructor default — 0,
    # which behaves as a non-sticky/non-sediment material).

    # Rectangular xy lattice sized so dx_m ≈ DX_TOPO (one marker per topo
    # column) and dy_m ≤ compaction cell height (one marker per cell, often
    # several). See MARKER_NX / MARKER_NY definitions at the top of this file.
    nx_m = MARKER_NX
    ny_m = MARKER_NY
    nm = nx_m * ny_m
    @assert nm <= nmarkers "Need $nm marker slots; only $nmarkers allocated. Bump nmarkers_cell_x/y."
    dx_m = XSIZE / Float64(nx_m)
    dy_m = YSIZE / Float64(ny_m)
    println(">> Marker grid: $nx_m × $ny_m (dx=$(round(dx_m, digits=1)) m, dy=$(round(dy_m, digits=1)) m)")
    println(">>   using $nm of $nmarkers slots")

    # Initialize all markers to a "neutral" matid 0 (basement-like). Then
    # overwrite the sticky-water region and the sediment basin region.
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
        if y < Y_SEALEVEL + STICKY_THICKNESS
            marker_matid[k] = matid_sticky_water
            marker_porosity_initial[k] = Float32(0.0)
            marker_decay_depth[k] = Float32(2_000.0)
            marker_max_burial_depth[k] = Float32(0.0)
        elseif x >= X_BASIN_MIN && x <= X_BASIN_MAX && y <= BASEMENT_DEPTH
            marker_matid[k] = matid_sediment
            marker_porosity_initial[k] = Float32(SEDIMENT_POROSITY_INITIAL)
            marker_decay_depth[k] = Float32(SEDIMENT_DECAY_DEPTH)
            # SEDIMENT_PRECOMPACTED selects whether the marker's prior burial
            # counts as already-compacted. See compact_sediment_mesh: cells
            # only compact when their bottom exceeds the average max-burial
            # depth of markers in that cell.
            marker_max_burial_depth[k] = SEDIMENT_PRECOMPACTED ?
                Float32(y - (Y_SEALEVEL + STICKY_THICKNESS)) : Float32(0.0)
        else
            # Basement / mantle below the sediment, or off-basin: leave matid=0.
            marker_matid[k] = Int16(0)
            marker_porosity_initial[k] = Float32(0.0)
            marker_decay_depth[k] = Float32(2_000.0)
            marker_max_burial_depth[k] = Float32(0.0)
        end
    end
    return nothing
end

function configure_drainage_basin!(model::ModelData)
    model.melting.parameters.extraction.ndrainage_basin.value = 1

    arrays = model.melting.arrays.extraction
    arrays.xstart_drainage.array[1]                  = 0.0
    arrays.xend_drainage.array[1]                    = XSIZE
    arrays.extrusion_volumes.array[1]                = TOTAL_EXTRUSION_VOLUME
    arrays.xmid_molten_zones.array[1]                = ERUPTION_X_MIN + ERUPTION_DOMAIN_WIDTH / 2.0
    arrays.width_molten_zones.array[1]               = ERUPTION_DOMAIN_WIDTH
    arrays.avg_shallow_partial_melt_xcoors.array[1]  = ERUPTION_X_MIN + ERUPTION_DOMAIN_WIDTH / 2.0
    arrays.avg_shallow_partial_melt_ycoors.array[1]  = Y_SEALEVEL + STICKY_THICKNESS
    return nothing
end

function build_sediment_thickness_initial(model::ModelData)::Vector{Float64}
    sediment_thickness = zeros(Float64, TOPONUM)
    for i in 1:TOPONUM
        x = (i - 1) * DX_TOPO
        if x >= X_BASIN_MIN && x <= X_BASIN_MAX
            sediment_thickness[i] = SEDIMENT_THICKNESS
        end
    end
    return sediment_thickness
end

# ── Snapshotting & plotting ─────────────────────────────────────────────────

struct State
    topo_gridx::Vector{Float64}
    topo_gridy::Vector{Float64}
    extrusion_thickness::Vector{Float64}
    marker_x::Vector{Float64}
    marker_y::Vector{Float64}
    marker_matid::Vector{Int16}
end

function snapshot_state(model::ModelData)::State
    gridt = model.topography.arrays.gridt.array
    return State(
        copy(gridt[1, :]),
        copy(gridt[2, :]),
        copy(gridt[7, :]),
        copy(model.markers.arrays.location.marker_x.array),
        copy(model.markers.arrays.location.marker_y.array),
        copy(model.markers.arrays.material.marker_matid.array),
    )
end

function print_diagnostics(pre::State, post::State)::Nothing
    lava = post.extrusion_thickness
    topo_delta = post.topo_gridy .- pre.topo_gridy
    marker_dy = post.marker_y .- pre.marker_y
    moved_markers = count(!=(0.0), marker_dy)

    println(">> Diagnostics:")
    println("   max lava thickness:        $(round(maximum(lava); digits=1)) m")
    println("   total lava deposit (∫dx):  $(round(sum(lava) * (pre.topo_gridx[2] - pre.topo_gridx[1]); digits=0)) m^2")
    println("   topo y change min/max:     $(round(minimum(topo_delta); digits=2)) / $(round(maximum(topo_delta); digits=2)) m")
    println("   markers displaced:         $moved_markers / $(length(marker_dy))")
    if moved_markers > 0
        nz = marker_dy[marker_dy .!= 0.0]
        println("   marker dy min/max:         $(round(minimum(nz); digits=2)) / $(round(maximum(nz); digits=2)) m")
    end
    return nothing
end

function plot_topography(pre::State, post::State)::Nothing
    fig = Figure(size = (1500, 500))
    ax = Axis(
        fig[1, 1];
        xlabel = "x (km)",
        ylabel = "y (m, depth-positive)",
        title = "Topography: pre vs post run_lava_flow_model",
    )
    lines!(ax, pre.topo_gridx ./ 1_000.0, pre.topo_gridy;
           color = :blue, linewidth = 2.0, label = "Topography (pre)")
    lines!(ax, post.topo_gridx ./ 1_000.0, post.topo_gridy;
           color = :red, linewidth = 2.0, linestyle = :dash, label = "Topography (post)")
    hlines!(ax, [Y_SEALEVEL]; color = :cyan, linestyle = :dot, label = "Sea level")
    xlims!(ax, 0.0, XSIZE / 1_000.0)
    ax.yreversed = true
    axislegend(ax; position = :rb)
    save("compaction_thick_sediment_topography.png", fig)
    return nothing
end

function plot_lava_thickness(post::State)::Nothing
    fig = Figure(size = (1500, 400))
    ax = Axis(
        fig[1, 1];
        xlabel = "x (km)",
        ylabel = "lava thickness (m)",
        title = "Total lava thickness deposited this timestep",
    )
    lines!(ax, post.topo_gridx ./ 1_000.0, post.extrusion_thickness;
           color = :darkorange, linewidth = 2.0,
           label = "Total lava thickness (gridt[7,:])")
    xlims!(ax, 0.0, XSIZE / 1_000.0)
    axislegend(ax; position = :rt)
    save("compaction_thick_sediment_lava.png", fig)
    return nothing
end

function plot_markers(model::ModelData, pre::State, post::State)::Nothing
    matid_types = model.materials.dicts.matid_types
    matid_sediment = matid_types["Sediment"][1]

    fig = Figure(size = (1500, 800))
    ax = Axis(
        fig[1, 1];
        xlabel = "x (km)",
        ylabel = "y (km, depth-positive)",
        title = "Sediment markers: pre vs post compaction-corrected extrusion",
    )

    pre_sed = pre.marker_matid .== matid_sediment
    post_sed = post.marker_matid .== matid_sediment

    scatter!(ax,
        pre.marker_x[pre_sed] ./ 1_000.0, pre.marker_y[pre_sed] ./ 1_000.0;
        markersize = 2.0, color = :blue, strokewidth = 0,
        label = "Sediment markers (pre)",
    )
    scatter!(ax,
        post.marker_x[post_sed] ./ 1_000.0, post.marker_y[post_sed] ./ 1_000.0;
        markersize = 2.0, color = :red, strokewidth = 0,
        label = "Sediment markers (post)",
    )
    lines!(ax,
        post.topo_gridx ./ 1_000.0, post.topo_gridy ./ 1_000.0;
        color = :black, linewidth = 2.0, label = "Topography (post)",
    )
    hlines!(ax, [Y_SEALEVEL ./ 1_000.0]; color = :cyan, linestyle = :dot, label = "Sea level")
    xlims!(ax, 0.0, XSIZE / 1_000.0)
    ylims!(ax, 0.0, YSIZE / 1_000.0)
    ax.yreversed = true
    axislegend(ax; position = :rb)
    save("compaction_thick_sediment_markers.png", fig)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    LavaFlowTestCompactionThickSediment.run_test()
end
