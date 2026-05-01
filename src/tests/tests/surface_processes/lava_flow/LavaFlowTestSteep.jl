module LavaFlowTestSteep

using CairoMakie
import EarthBox.SurfaceProcesses.LavaFlowManager.LavaFlowSolverManager.MakeFlow: make_flow

"""
Demonstrate the steep-wall basin-trapping bug using topography that replicates
the SDR experiment geometry: a broad, gently-sloped rift valley with a narrow,
near-vertical-walled axial slot at the centre (as visible in SDR case 10).

Two cases with identical eruption volume and multiple pulses (use_single_pulse=false):

  1. Composite rift (broad + slot): lava erupts into the narrow slot and is
     trapped by the steep slot walls — it cannot spread to the broad rift floor.
     Overflow requires lava > slot_depth_extra (~1700 m) at the slot edge.

  2. Broad rift only (no slot): the same volume easily overflows the gentle broad
     rift walls (threshold ~5 m) and spreads across the flanks.
"""
function run_test()
    xsize = 300_000.0
    xnum  = 601
    topo_gridx = collect(range(0, xsize, length=xnum))
    dx = topo_gridx[2] - topo_gridx[1]   # 500 m

    x_eruption         = xsize / 2.0
    residual_thickness = 30.0
    tolerance          = 1e-4
    nmax               = 1000
    decimation_factor  = 1

    # Broad rift parameters (outer structure)
    broad_half_width          = 30_000.0   # m  (±30 km, matching SDR rift scale)
    broad_depth_at_slot_edge  = 300.0      # m  (rift floor depth just outside slot)

    # Axial slot parameters (inner narrow trough)
    slot_half_width  = 2_000.0             # m  (±2 km flat floor, 8 nodes)
    slot_floor_depth = 2_000.0             # m  (total depth of slot floor)
    slot_wall_width  = dx                  # m  (1 node → wall ≈ 73°)
    slot_depth_extra = slot_floor_depth - broad_depth_at_slot_edge  # 1700 m

    wall_angle_deg = atand(slot_depth_extra / slot_wall_width)

    # Overflow thresholds
    broad_rift_step  = broad_depth_at_slot_edge /
                       ((broad_half_width - slot_half_width - slot_wall_width) / dx)
    slot_threshold   = slot_depth_extra    # lava must exceed this at slot edge

    # Eruption volume: ~15 % of slot fill capacity
    # Slot fill capacity = slot_depth_extra × (slot_half_width × 2) = 6.8M m²
    total_flow_volume = 1_000_000.0        # m²

    println("=== LavaFlowTestSteep ===")
    println("Broad rift half-width      : ", broad_half_width/1000, " km")
    println("Broad rift depth at edge   : ", broad_depth_at_slot_edge, " m")
    println("Broad rift step per node   : ≈ ", round(broad_rift_step, digits=1), " m  (easy overflow)")
    println("Slot half-width            : ", slot_half_width/1000, " km  (",
            round(Int, slot_half_width/dx), " floor nodes)")
    println("Slot floor depth           : ", slot_floor_depth, " m")
    println("Slot extra depth           : ", slot_depth_extra, " m")
    println("Slot wall angle            : ≈ ", round(wall_angle_deg, digits=0), "°")
    println("Slot overflow threshold    : > ", slot_threshold, " m at slot edge")
    println("Total flow volume          : ", total_flow_volume, " m²")
    println("Slot fill capacity         : ", slot_depth_extra * slot_half_width * 2, " m²")
    println("Flow / slot-fill ratio     : ",
            round(total_flow_volume / (slot_depth_extra * slot_half_width * 2) * 100, digits=1), "%")
    println()

    # ── Case 1: composite rift (broad gentle rift + narrow steep axial slot) ────
    topo_composite = make_topoy_composite_rift(
        topo_gridx, x_eruption,
        slot_half_width, slot_floor_depth, slot_wall_width,
        broad_depth_at_slot_edge, broad_half_width
    )
    lava_composite = zeros(xnum)
    make_flow(topo_gridx, topo_composite, lava_composite, total_flow_volume,
              residual_thickness, x_eruption;
              decimation_factor=decimation_factor, tolerance=tolerance, nmax=nmax,
              use_single_pulse=false)

    print_diagnostics("Composite rift (broad + steep slot)",
                      lava_composite, x_eruption, slot_half_width, total_flow_volume, dx)

    # ── Case 2: broad rift only, no slot ────────────────────────────────────────
    topo_broad = make_topoy_broad_rift(topo_gridx, x_eruption,
                                       broad_half_width, broad_depth_at_slot_edge)
    lava_broad = zeros(xnum)
    make_flow(topo_gridx, topo_broad, lava_broad, total_flow_volume,
              residual_thickness, x_eruption;
              decimation_factor=decimation_factor, tolerance=tolerance, nmax=nmax,
              use_single_pulse=false)

    print_diagnostics("Broad rift only (no slot)",
                      lava_broad, x_eruption, slot_half_width, total_flow_volume, dx)

    plot_surface_comparison(topo_gridx, topo_composite, lava_composite,
                            topo_broad, lava_broad,
                            slot_threshold, broad_rift_step)
    plot_thickness_comparison(topo_gridx, lava_composite, lava_broad)
end

# ── Topography builders ──────────────────────────────────────────────────────

"""
Composite SDR-style topography: broad gentle rift with a narrow steep axial slot.

  dist ∈ [0, slot_half_width]                  → flat slot floor at slot_floor_depth
  dist ∈ (slot_half_width, +slot_wall_width]   → steep wall: slot_floor_depth → broad_depth_at_slot_edge
  dist ∈ (slot+wall, broad_half_width]         → gentle broad rift floor tapering to 0
  dist > broad_half_width                      → flat flanks at 0
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
Broad rift only — same outer structure but no axial slot.
Linear V-shape from broad_depth_at_center to 0 over broad_half_width.
Per-node elevation step is small so lava overflows easily.
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

# ── Diagnostics ─────────────────────────────────────────────────────────────

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
    total_placed  = sum(lava_thickness) * dx

    n = length(lava_thickness)
    center_node = eruption_node
    slot_nodes  = floor(Int, slot_half_width / dx)
    i_lo = max(1, center_node - slot_nodes)
    i_hi = min(n, center_node + slot_nodes)

    in_slot_volume  = sum(lava_thickness[i_lo:i_hi]) * dx
    out_slot_volume = total_placed - in_slot_volume
    out_frac        = out_slot_volume / max(total_placed, 1e-10) * 100

    spread_width_m = count(t -> t > 0.1, lava_thickness) * dx

    println("── $label ──")
    println("  Max lava thickness          : ", round(max_thickness, digits=1), " m")
    println("  Thickness at eruption node  : ", round(lava_at_source, digits=1), " m")
    println("  Volume inside slot          : ", round(in_slot_volume, digits=0), " m²")
    println("  Volume outside slot         : ", round(out_slot_volume, digits=0),
            " m²  (", round(out_frac, digits=1), "% of total)")
    println("  Lateral spread width        : ", round(spread_width_m / 1000, digits=1), " km")
    println("  Total volume placed         : ", round(total_placed, digits=0),
            " m²  (input: ", round(Int, total_flow_volume), " m²)")
    println()
    return nothing
end

# ── Plots ────────────────────────────────────────────────────────────────────

function plot_surface_comparison(
    topo_gridx, topo_composite, lava_composite, topo_broad, lava_broad,
    slot_threshold, broad_rift_step
)::Nothing
    xkm = topo_gridx ./ 1000.0

    fig = Figure(size = (1200, 700))
    ax1 = Axis(
        fig[1, 1];
        xlabel = "X (km)", ylabel = "Depth (m)",
        title = "Composite rift — lava TRAPPED in slot (threshold $(round(Int,slot_threshold)) m)",
    )
    lines!(ax1, xkm, topo_composite; color = :black, linewidth = 2, label = "Topography")
    lines!(ax1, xkm, topo_composite .- lava_composite; color = :red, linewidth = 2, label = "Lava surface")
    ax1.yreversed = true
    axislegend(ax1)

    ax2 = Axis(
        fig[2, 1];
        xlabel = "X (km)", ylabel = "Depth (m)",
        title = "Broad rift only — lava ESCAPES (threshold ≈ $(round(broad_rift_step, digits=1)) m)",
    )
    lines!(ax2, xkm, topo_broad; color = :black, linewidth = 2, label = "Topography")
    lines!(ax2, xkm, topo_broad .- lava_broad; color = :blue, linewidth = 2, label = "Lava surface")
    ax2.yreversed = true
    axislegend(ax2)

    save("lava_flow_steep_vs_gentle_surface.png", fig)

    # Zoomed view centred on rift axis (±50 km)
    zoom_km = 50.0
    x_centre_km = (topo_gridx[end] / 2.0) / 1000.0
    xlims_zoom = (x_centre_km - zoom_km, x_centre_km + zoom_km)

    fig_z = Figure(size = (1200, 700))
    axz1 = Axis(
        fig_z[1, 1];
        xlabel = "X (km)", ylabel = "Depth (m)",
        title = "Composite rift — zoomed (±$(round(Int,zoom_km)) km)",
    )
    lines!(axz1, xkm, topo_composite; color = :black, linewidth = 2, label = "Topography")
    lines!(axz1, xkm, topo_composite .- lava_composite; color = :red, linewidth = 2, label = "Lava surface")
    xlims!(axz1, xlims_zoom...)
    axz1.yreversed = true
    axislegend(axz1)

    axz2 = Axis(
        fig_z[2, 1];
        xlabel = "X (km)", ylabel = "Depth (m)",
        title = "Broad rift only — zoomed (±$(round(Int,zoom_km)) km)",
    )
    lines!(axz2, xkm, topo_broad; color = :black, linewidth = 2, label = "Topography")
    lines!(axz2, xkm, topo_broad .- lava_broad; color = :blue, linewidth = 2, label = "Lava surface")
    xlims!(axz2, xlims_zoom...)
    axz2.yreversed = true
    axislegend(axz2)

    save("lava_flow_steep_vs_gentle_surface_zoom.png", fig_z)
    return nothing
end

function plot_thickness_comparison(
    topo_gridx, lava_composite, lava_broad
)::Nothing
    xkm = topo_gridx ./ 1000.0
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

    save("lava_flow_steep_vs_gentle_thickness.png", fig)
    return nothing
end

end # module
