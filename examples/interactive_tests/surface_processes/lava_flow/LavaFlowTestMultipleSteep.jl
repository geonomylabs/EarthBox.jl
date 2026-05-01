module LavaFlowTestMultipleSteep

using CairoMakie
import EarthBox.ConversionFuncs: mm_per_yr_to_meters_per_seconds as mm_yr_to_m_s
import EarthBox.ConversionFuncs: years_to_seconds
import EarthBox.ConversionFuncs: meters_per_year_to_meters_per_seconds as m_yr_to_m_s
import EarthBox.ConversionFuncs: meters_squared_per_year_to_meters_squared_per_second as m2_yr_to_m2_s
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.SurfaceProcesses.LavaFlowManager: calculate_number_of_flows
import EarthBox.SurfaceProcesses.LavaFlowManager: calculate_characteristic_volume_per_flow
import EarthBox.SurfaceProcesses.LavaFlowManager.LavaFlowSolverManager: LavaFlowSolver
import EarthBox.SurfaceProcesses.LavaFlowManager.LavaFlowSolverManager: extrude_magma

"""
Test that replicates the actual SDR model's lava flow pipeline (LavaFlowSolver +
extrude_magma) on two contrasting rift geometries:

  1. Composite rift: broad gentle rift with a narrow steep axial slot (SDR-like).
     Multiple flows erupt into the slot and are trapped by the steep walls,
     producing a tall narrow edifice.

  2. Broad rift only: same outer geometry but no axial slot.
     The same total volume spreads widely across the rift floor.

Both cases use submarine residual thickness (50 m) and deterministic eruption at
the rift axis centre, matching early SDR conditions.
"""
function run_test()
    xsize = 500_000.0
    dx    = 100.0
    xnum  = floor(Int, xsize / dx) + 1
    topo_gridx = collect(range(0.0, xsize, length=xnum))

    x_eruption = xsize / 2.0   # rift axis at domain centre

    # Broad rift geometry
    broad_half_width         = 30_000.0   # m
    broad_depth_at_slot_edge = 1000.0      # m

    # Axial slot geometry
    slot_half_width  = 2_000.0            # m  (±2 km, 4 floor nodes at dx=500m)
    slot_floor_depth = 2_000.0            # m
    slot_wall_width  = dx                 # m  (1 node → ~74° walls)
    slot_depth_extra = slot_floor_depth - broad_depth_at_slot_edge

    # Eruption domain: centred on rift axis, spanning the slot floor
    eruption_location_x_min = x_eruption - slot_half_width
    width_eruption_domain   = 2_500.0 #slot_half_width * 2.0

    # Eruption parameters
    total_extrusion_volume             = 16_000_000.0   # m²
    characteristic_flow_length         = 60_000.0
    residual_lava_thickness_subaerial  = 30.0           # m
    residual_lava_thickness_submarine  = 30.0           # m # 50 m
    characteristic_volume_per_flow = characteristic_flow_length * residual_lava_thickness_submarine
    number_of_flows = max(1, floor(Int, total_extrusion_volume / characteristic_volume_per_flow))
    println("characteristic_volume_per_flow: ", characteristic_volume_per_flow)
    println("total_extrusion_volume: ", total_extrusion_volume)
    println("number_of_flows: ", number_of_flows)
    y_sealevel                         = 0.0            # everything submarine
    decimation_factor                  = 8

    # Dummy decompaction parameters (compaction correction disabled)
    decompaction_params = SedimentTransportParameters(
        m2_yr_to_m2_s(0.25),
        m_yr_to_m_s(1.0),
        1e-3,
        m2_yr_to_m2_s(1e2),
        2000.0,
        years_to_seconds(50_000.0),
        years_to_seconds(5.0e6),
        0.4,
        1 / 2000.0
    )

    sediment_initial = zeros(xnum)

    broad_rift_step = broad_depth_at_slot_edge /
                      ((broad_half_width - slot_half_width - slot_wall_width) / dx)

    println("=== LavaFlowTestMultipleSteep ===")
    println("Grid: $(xnum) nodes, dx=$(dx) m, domain=$(xsize/1000) km")
    println("Slot half-width        : $(slot_half_width/1000) km")
    println("Slot floor depth       : $(slot_floor_depth) m")
    println("Slot overflow threshold: $(slot_depth_extra) m")
    println("Broad rift step/node   : ≈ $(round(broad_rift_step, digits=1)) m")
    println("Number of flows        : $(number_of_flows)")
    println("Volume per flow        : $(total_extrusion_volume / number_of_flows) m²")
    println("Submarine residual     : $(residual_lava_thickness_submarine) m")
    println()

    # ── Case 1: composite rift ───────────────────────────────────────────────
    topo_composite = make_topoy_composite_rift(
        topo_gridx, x_eruption,
        slot_half_width, slot_floor_depth, slot_wall_width,
        broad_depth_at_slot_edge, broad_half_width
    )
    topo_composite_initial = copy(topo_composite)

    solver_composite = LavaFlowSolver(
        topo_gridx, topo_composite, sediment_initial,
        eruption_location_x_min, width_eruption_domain,
        total_extrusion_volume, number_of_flows,
        residual_lava_thickness_subaerial, residual_lava_thickness_submarine,
        y_sealevel, decompaction_params;
        use_random_eruption_location=true,
        use_compaction_correction=false,
        decimation_factor=decimation_factor
    )

    extrude_magma(solver_composite)

    print_diagnostics("Composite rift (broad + steep slot)",
                      solver_composite.total_lava_thickness,
                      x_eruption, slot_half_width, total_extrusion_volume, dx)

    # ── Case 2: broad rift only ──────────────────────────────────────────────
    topo_broad = make_topoy_broad_rift(
        topo_gridx, x_eruption, broad_half_width, broad_depth_at_slot_edge
    )
    topo_broad_initial = copy(topo_broad)

    solver_broad = LavaFlowSolver(
        topo_gridx, topo_broad, sediment_initial,
        eruption_location_x_min, width_eruption_domain,
        total_extrusion_volume, number_of_flows,
        residual_lava_thickness_subaerial, residual_lava_thickness_submarine,
        y_sealevel, decompaction_params;
        use_random_eruption_location=true,
        use_normal_eruption_location=true,
        use_compaction_correction=false,
        decimation_factor=decimation_factor
    )

    extrude_magma(solver_broad)

    print_diagnostics("Broad rift only (no slot)",
                      solver_broad.total_lava_thickness,
                      x_eruption, slot_half_width, total_extrusion_volume, dx)

    plot_surface_comparison(
        topo_gridx,
        topo_composite_initial, solver_composite.total_lava_thickness,
        topo_broad_initial,     solver_broad.total_lava_thickness,
        slot_depth_extra, broad_rift_step
    )
    plot_thickness_comparison(
        topo_gridx,
        solver_composite.total_lava_thickness,
        solver_broad.total_lava_thickness
    )
end

# ── Topography builders ──────────────────────────────────────────────────────

"""
Composite SDR-style topography: broad gentle rift with a narrow steep axial slot.

  dist ∈ [0, slot_half_width]                → flat slot floor at slot_floor_depth
  dist ∈ (slot_half_width, +slot_wall_width] → steep wall: slot_floor_depth → broad_depth_at_slot_edge
  dist ∈ (slot+wall, broad_half_width]       → gentle broad rift tapering to 0
  dist > broad_half_width                    → flat flanks at 0
"""
function make_topoy_composite_rift(
    topo_gridx::Vector{Float64},
    valley_center::Float64,
    slot_half_width::Float64,
    slot_floor_depth::Float64,
    slot_wall_width::Float64,
    broad_depth_at_slot_edge::Float64,
    broad_half_width::Float64
)::Vector{Float64}
    topo_gridy = zeros(length(topo_gridx))
    broad_outer_width = broad_half_width - slot_half_width - slot_wall_width
    for i in eachindex(topo_gridx)
        dist = abs(topo_gridx[i] - valley_center)
        if dist <= slot_half_width
            topo_gridy[i] = slot_floor_depth
        elseif dist <= slot_half_width + slot_wall_width
            frac = (dist - slot_half_width) / slot_wall_width
            topo_gridy[i] = slot_floor_depth * (1 - frac) + broad_depth_at_slot_edge * frac
        elseif dist <= broad_half_width
            frac = (dist - slot_half_width - slot_wall_width) / broad_outer_width
            topo_gridy[i] = broad_depth_at_slot_edge * (1 - frac)
        end
    end
    return topo_gridy
end

"""
Broad rift only — linear V from broad_depth_at_center at the axis to 0 at broad_half_width.
Per-node step is small so lava overflows easily.
"""
function make_topoy_broad_rift(
    topo_gridx::Vector{Float64},
    valley_center::Float64,
    broad_half_width::Float64,
    broad_depth_at_center::Float64
)::Vector{Float64}
    topo_gridy = zeros(length(topo_gridx))
    for i in eachindex(topo_gridx)
        dist = abs(topo_gridx[i] - valley_center)
        if dist <= broad_half_width
            topo_gridy[i] = broad_depth_at_center * (1.0 - dist / broad_half_width)
        end
    end
    return topo_gridy
end

# ── Diagnostics ──────────────────────────────────────────────────────────────

function print_diagnostics(
    label::String,
    lava_thickness::Vector{Float64},
    x_eruption::Float64,
    slot_half_width::Float64,
    total_flow_volume::Float64,
    dx::Float64
)::Nothing
    eruption_node = floor(Int, x_eruption / dx) + 1
    max_thickness = maximum(lava_thickness)
    lava_at_source = lava_thickness[eruption_node]
    total_placed   = sum(lava_thickness) * dx

    n = length(lava_thickness)
    slot_nodes = floor(Int, slot_half_width / dx)
    i_lo = max(1, eruption_node - slot_nodes)
    i_hi = min(n, eruption_node + slot_nodes)

    in_slot_volume  = sum(lava_thickness[i_lo:i_hi]) * dx
    out_slot_volume = total_placed - in_slot_volume
    out_frac        = out_slot_volume / max(total_placed, 1e-10) * 100

    spread_width_m = count(t -> t > 0.1, lava_thickness) * dx

    println("── $label ──")
    println("  Max lava thickness         : $(round(max_thickness, digits=1)) m")
    println("  Thickness at eruption node : $(round(lava_at_source, digits=1)) m")
    println("  Volume inside slot         : $(round(in_slot_volume, digits=0)) m²")
    println("  Volume outside slot        : $(round(out_slot_volume, digits=0)) m²  ($(round(out_frac, digits=1))% of total)")
    println("  Lateral spread width       : $(round(spread_width_m / 1000, digits=1)) km")
    println("  Total volume placed        : $(round(total_placed, digits=0)) m²  (input: $(round(Int, total_flow_volume)) m²)")
    println()
    return nothing
end

# ── Plots ─────────────────────────────────────────────────────────────────────

function plot_surface_comparison(
    topo_gridx,
    topo_composite, lava_composite,
    topo_broad, lava_broad,
    slot_threshold, broad_rift_step
)::Nothing
    xkm = topo_gridx ./ 1000.0

    fig = Figure(size = (1200, 700))
    ax1 = Axis(
        fig[1, 1];
        xlabel = "X (km)", ylabel = "Depth (m)",
        title = "Composite rift — lava TRAPPED in slot (threshold $(round(Int,slot_threshold)) m)",
    )
    lines!(ax1, xkm, topo_composite; color = :black, linewidth = 2, label = "Initial topo")
    lines!(ax1, xkm, topo_composite .- lava_composite; color = :red, linewidth = 2, label = "Lava surface")
    ax1.yreversed = true
    axislegend(ax1)

    ax2 = Axis(
        fig[2, 1];
        xlabel = "X (km)", ylabel = "Depth (m)",
        title = "Broad rift only — lava ESCAPES (threshold ≈ $(round(broad_rift_step, digits=1)) m/node)",
    )
    lines!(ax2, xkm, topo_broad; color = :black, linewidth = 2, label = "Initial topo")
    lines!(ax2, xkm, topo_broad .- lava_broad; color = :blue, linewidth = 2, label = "Lava surface")
    ax2.yreversed = true
    axislegend(ax2)

    save("lava_flow_multiple_steep_surface.png", fig)

    # Zoomed view ±50 km around rift axis
    zoom_km     = 50.0
    x_centre_km = (topo_gridx[end] / 2.0) / 1000.0
    xlims_zoom  = (x_centre_km - zoom_km, x_centre_km + zoom_km)

    fig_z = Figure(size = (1200, 700))
    axz1 = Axis(
        fig_z[1, 1];
        xlabel = "X (km)", ylabel = "Depth (m)",
        title = "Composite rift — zoomed (±$(round(Int,zoom_km)) km)",
    )
    lines!(axz1, xkm, topo_composite; color = :black, linewidth = 2, label = "Initial topo")
    lines!(axz1, xkm, topo_composite .- lava_composite; color = :red, linewidth = 2, label = "Lava surface")
    xlims!(axz1, xlims_zoom...)
    axz1.yreversed = true
    axislegend(axz1)

    axz2 = Axis(
        fig_z[2, 1];
        xlabel = "X (km)", ylabel = "Depth (m)",
        title = "Broad rift only — zoomed (±$(round(Int,zoom_km)) km)",
    )
    lines!(axz2, xkm, topo_broad; color = :black, linewidth = 2, label = "Initial topo")
    lines!(axz2, xkm, topo_broad .- lava_broad; color = :blue, linewidth = 2, label = "Lava surface")
    xlims!(axz2, xlims_zoom...)
    axz2.yreversed = true
    axislegend(axz2)

    save("lava_flow_multiple_steep_surface_zoom.png", fig_z)
    return nothing
end

function plot_thickness_comparison(
    topo_gridx, lava_composite, lava_broad
)::Nothing
    xkm       = topo_gridx ./ 1000.0
    bar_width = xkm[2] - xkm[1]

    fig = Figure(size = (1200, 700))
    ax1 = Axis(
        fig[1, 1];
        xlabel = "X (km)", ylabel = "Thickness (m)",
        title = "Composite rift — lava thickness (edifice in slot)",
    )
    barplot!(ax1, xkm, lava_composite;
             width = bar_width, color = :red, strokewidth = 0, label = "Lava thickness")
    axislegend(ax1)

    ax2 = Axis(
        fig[2, 1];
        xlabel = "X (km)", ylabel = "Thickness (m)",
        title = "Broad rift only — lava thickness (spread)",
    )
    barplot!(ax2, xkm, lava_broad;
             width = bar_width, color = :blue, strokewidth = 0, label = "Lava thickness")
    axislegend(ax2)

    save("lava_flow_multiple_steep_thickness.png", fig)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    LavaFlowTestMultipleSteep.run_test()
end
