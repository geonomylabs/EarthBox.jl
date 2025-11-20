module LavaFlowTestSingle

import Plots
import EarthBox.SurfaceProcesses.LavaFlowManager.LavaFlowSolverManager.MakeFlow: make_flow
import EarthBox.SurfaceProcesses.LavaFlowManager.LavaFlowSolverManager.MakeFlow.LavaFlowPulse: radiate_indices

function test_radiate_indices()
    sorted_indices = radiate_indices(11, 5)
    println("sorted_indices: ", sorted_indices)
end

function run_test()
    xnum = 5001 # Number of grid points in the x direction for topography grid
    xsize = 500_000.0 # Size of the topography grid in the x direction (meters)
    flow_grid_decimation_factor = 10
    characteristic_flow_length = 100_000.0 # meters
    residual_lava_thickness = 30.0 # Residual lava thickness (meters)
    total_flow_volume = characteristic_flow_length * residual_lava_thickness * 3.0
    println(">> total_flow_volume: ", total_flow_volume)
    x_eruption_location = 250_000.0 # Eruption location in the x direction (meters)
    topography_dip = 0.1 # Dip of the topography in degrees
    basin_start = 300_000.0 # meters
    basin_width = 50_000.0 # meters
    basin_depth = 50.0 # meters
    tolerance = 1e-4
    nmax = 1000
    use_single_pulse = false

    topo_gridx = collect(range(0, xsize, length=xnum))
    lava_thickness = zeros(xnum)
    topo_gridy = make_topoy_ridge(topo_gridx, topography_dip)
    add_basin(topo_gridx, topo_gridy, basin_start, basin_width, basin_depth)

    t1 = time()
    make_flow(
        topo_gridx, topo_gridy, lava_thickness, total_flow_volume,
        residual_lava_thickness, x_eruption_location;
        decimation_factor=flow_grid_decimation_factor,
        tolerance=tolerance, nmax=nmax, use_single_pulse=use_single_pulse
    )
    t2 = time()
    println("Time to run make_flow (seconds): ", t2 - t1)

    plot_lava_thickness(99999, 99999, 99999, topo_gridx, lava_thickness)
    plot_surface_and_topography(
        99999, 99999, 99999, topo_gridx, topo_gridy, lava_thickness;
        figsize=(5, 2), extension=".png"
    )
end

function make_topoy_slope(
    topo_gridx::Vector{Float64},
    dip_degrees::Float64
)::Vector{Float64}
    xnum = length(topo_gridx)
    topo_gridy = zeros(xnum)
    for i in 1:xnum
        topo_gridy[i] = topo_gridx[i] * tan(deg2rad(dip_degrees))
    end
    return topo_gridy
end

function make_topoy_ridge(
    topo_gridx::Vector{Float64},
    dip_degrees::Float64
)::Vector{Float64}
    xnum = length(topo_gridx)
    topo_gridy = zeros(xnum)
    xnum_mid = floor(Int, xnum/2)
    xmid = topo_gridx[xnum_mid]
    for i in 1:xnum_mid
        topo_gridy[i] = (xmid - topo_gridx[i]) * tan(deg2rad(dip_degrees))
    end
    for i in xnum_mid+1:xnum
        topo_gridy[i] = (topo_gridx[i] - xmid) * tan(deg2rad(dip_degrees))
    end
    return topo_gridy
end

function add_basin(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    basin_start::Float64,
    basin_width::Float64,
    basin_amp::Float64
)::Nothing
    dx_smooth = basin_width/5.0
    for i in 1:length(topo_gridx)
        x = topo_gridx[i]
        if x > basin_start && x < basin_start + basin_width
            if basin_start + dx_smooth <= x <= basin_start + basin_width - dx_smooth
                topo_gridy[i] = topo_gridy[i] + basin_amp
            elseif x < basin_start + dx_smooth
                topo_gridy[i] = topo_gridy[i] + basin_amp * (x - basin_start)/dx_smooth
            elseif x > basin_start + basin_width - dx_smooth
                topo_gridy[i] = topo_gridy[i] + basin_amp * 
                    (basin_start + basin_width - x)/dx_smooth
            end
        end
    end
    return nothing
end

function print_iteration_information(
    i::Int,
    topo_left::Float64,
    topo_active::Float64,
    topo_right::Float64,
    y_surface_left::Float64,
    y_surface_active::Float64,
    y_surface_right::Float64,
    delta_elevation_left::Float64,
    delta_elevation_right::Float64,
    delta_elevation_total::Float64,
    max_potential_outflow_thickness::Float64,
    total_potential_outflow_thickness::Float64,
    remaining_excess_thickness::Float64,
    outflow_thickness_to_left::Float64,
    outflow_thickness_to_right::Float64,
    old_thickness::Float64,
    lava_thickness::Vector{Float64}
)::Nothing
    println("i: ", i)
    println("topo_left: ", topo_left)
    println("topo_active: ", topo_active)
    println("topo_right: ", topo_right)
    println("y_surface_left: ", y_surface_left)
    println("y_surface_active: ", y_surface_active)
    println("y_surface_right: ", y_surface_right)
    println("delta_elevation_left: ", delta_elevation_left)
    println("delta_elevation_right: ", delta_elevation_right)
    println("delta_elevation_total: ", delta_elevation_total)
    println("max_potential_outflow_thickness: ", max_potential_outflow_thickness)
    println("total_potential_outflow_thickness: ", total_potential_outflow_thickness)
    println("remaining_excess_thickness: ", remaining_excess_thickness)
    println("outflow_thickness_to_left: ", outflow_thickness_to_left)
    println("outflow_thickness_to_right: ", outflow_thickness_to_right)
    println("old lava_thickness: ", old_thickness)
    println("new lava_thickness[i]: ", lava_thickness[i])
    return nothing
end

function plot_surface_and_topography(
    pulse_id::Int,
    icount::Int,
    i::Int,
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    lava_thickness::Vector{Float64};
    figsize::Tuple{Int,Int}=(10, 2),
    extension::String=".png"
)::Nothing
    dpi = 150
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)
    p = Plots.plot(
        topo_gridx, topo_gridy .- lava_thickness,
        label="Surface",
        size=figsize_pixels,
        margin=5Plots.mm
    )
    Plots.plot!(p, topo_gridx, topo_gridy, label="Topography")
    #Plots.ylims!(p, (100.0, -100.0))
    Plots.title!("Pulse $pulse_id Iteration $icount Index $i")
    Plots.xlabel!("X(m)")
    Plots.ylabel!("Y(m)")
    Plots.yaxis!(:flip)
    Plots.savefig(p, "surface_p$(pulse_id)_icount$(icount)$(extension)")
    return nothing
end

function plot_lava_thickness(
    pulse_id::Int,
    icount::Int,
    i::Int,
    topo_gridx::Vector{Float64},
    lava_thickness::Vector{Float64}
)::Nothing
    bar_width = topo_gridx[2] - topo_gridx[1]
    p = Plots.bar(
        topo_gridx, lava_thickness,
        bar_width=bar_width,
        label="Lava Thickness",
        margin=5Plots.mm,
        linewidth=0,  # Turn off bar edges
        fillcolor=:blue  # Color all bars blue
    )
    Plots.title!("Lava Thickness at pulse $pulse_id Iteration $icount Index $i")
    Plots.xlabel!("Position")
    Plots.ylabel!("Lava Thickness")
    Plots.savefig(p, "lava_thickness_p$(pulse_id)_icount$(icount).png")
    return nothing
end

end # module 