module DebugPlots

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Arrays: ArrayUtils
import EarthBox.ConversionFuncs: get_factor_cm_yr_to_m_s
import CairoMakie

function plot_interpolated_temperature_tk1(model::ModelData, msg::String)::Nothing
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    ntimestep_str = string(ntimestep)
    gridx_b = model.grids.arrays.basic.gridx_b.array ./ 1000
    gridy_b = model.grids.arrays.basic.gridy_b.array ./ 1000
    tk1 = model.heat_equation.arrays.temperature.tk1.array

    ArrayUtils.print_min_max("tk1", tk1)

    dpi = 150
    figsize = (15, 5)
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)

    xspacing = 50
    yspacing = 25

    fig = CairoMakie.Figure(size = figsize_pixels)
    ax = CairoMakie.Axis(
        fig[1, 1];
        xlabel = "x (km)", ylabel = "y (km)",
        title = "Temperature (K)",
        titlesize = 25,
        xlabelsize = 25, ylabelsize = 25,
        xticklabelsize = 20, yticklabelsize = 20,
        xticks = minimum(gridx_b):xspacing:maximum(gridx_b),
        yticks = minimum(gridy_b):yspacing:maximum(gridy_b),
    )
    CairoMakie.xlims!(ax, minimum(gridx_b), maximum(gridx_b))
    CairoMakie.ylims!(ax, minimum(gridy_b), maximum(gridy_b))
    ax.yreversed = true
    hm = CairoMakie.heatmap!(
        ax, gridx_b, gridy_b, permutedims(tk1);
        colormap = :viridis,
    )
    CairoMakie.Colorbar(fig[1, 2], hm)
    CairoMakie.save("tk1_2d_grid_$(ntimestep_str).png", fig)

    return nothing
end

function plot_interpolated_basic_grid_velocity(model::ModelData, msg::String)::Nothing
    println("DEBUG: plotting interpolated basic grid velocity from ", msg)

    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    ntimestep_str = string(ntimestep)
    
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array

    # Vx-grid: y coordinates of vx grid nodes; note that x-coordinates are gridx_b
    #gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array 

    # Vy-grid: x coordinates of vy grid nodes; note that y coordinates are gridy_b
    #gridx_vy = model.grids.arrays.staggered_vy.gridx_vy.array 

    #vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    #vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array

    vxb = model.stokes_continuity.arrays.basic_grid_velocity.vxb.array
    vyb = model.stokes_continuity.arrays.basic_grid_velocity.vyb.array

    figsize = (30, 30)
    figsize_pixels = (figsize[1] * 100, figsize[2] * 100)

    cm_yr_to_m_s = get_factor_cm_yr_to_m_s()
    m_s_to_cm_yr = 1.0 / cm_yr_to_m_s

    # Convert velocities from m/s to cm/yr
    vxb_cm_yr = vxb .* m_s_to_cm_yr
    vyb_cm_yr = vyb .* m_s_to_cm_yr

    # Calculate scaling factor for quiver plot
    scale = 1.0 / maximum(sqrt.(vxb_cm_yr.^2 + vyb_cm_yr.^2))*1e4

    x = zeros(length(gridx_b)*length(gridy_b))
    y = zeros(length(gridx_b)*length(gridy_b))
    vx = zeros(length(gridx_b)*length(gridy_b))
    vy = zeros(length(gridx_b)*length(gridy_b))
    
    xnum = length(gridx_b)
    ynum = length(gridy_b)

    index = 1
    for j in 1:xnum
        for i in 1:ynum
            x[index] = gridx_b[j]
            y[index] = gridy_b[i]
            vx[index] = vxb_cm_yr[i, j] * scale
            vy[index] = vyb_cm_yr[i, j] * scale
            index += 1
        end
    end
    
    fig = CairoMakie.Figure(size = figsize_pixels)
    ax = CairoMakie.Axis(
        fig[1, 1];
        xlabel = "x (m)", ylabel = "y (m)",
        title = "Basic Grid Velocity Field (Normalized)",
        titlesize = 15,
        xlabelsize = 12, ylabelsize = 12,
        xticklabelsize = 10, yticklabelsize = 10,
        xticks = minimum(gridx_b):500:maximum(gridx_b),
        yticks = minimum(gridy_b):500:maximum(gridy_b),
    )
    CairoMakie.xlims!(ax, minimum(gridx_b), maximum(gridx_b))
    CairoMakie.ylims!(ax, minimum(gridy_b), maximum(gridy_b))
    CairoMakie.arrows!(ax, x, y, vx, vy; lengthscale = 1.0, color = :black)
    CairoMakie.save("basic_grid_velocity_field_$(ntimestep_str).png", fig)

    fig2 = CairoMakie.Figure(size = figsize_pixels)
    ax2 = CairoMakie.Axis(
        fig2[1, 1];
        xlabel = "x (m)", ylabel = "y (m)",
        title = "Vx (cm/yr x 1e4) Basic Grid",
        titlesize = 25,
        xlabelsize = 25, ylabelsize = 25,
        xticklabelsize = 20, yticklabelsize = 20,
        xticks = minimum(gridx_b):500:maximum(gridx_b),
        yticks = minimum(gridy_b):500:maximum(gridy_b),
    )
    CairoMakie.xlims!(ax2, minimum(gridx_b), maximum(gridx_b))
    CairoMakie.ylims!(ax2, minimum(gridy_b), maximum(gridy_b))
    hm2 = CairoMakie.heatmap!(
        ax2, gridx_b, gridy_b, permutedims(vxb_cm_yr) .* 1e4;
        colormap = :viridis,
    )
    CairoMakie.Colorbar(fig2[1, 2], hm2)
    plot_basic_grid_lines(ax2, model)
    plot_marker_overlay(ax2, model)
    CairoMakie.save("vx_basic_cm_yr_$(ntimestep_str).png", fig2)

    fig3 = CairoMakie.Figure(size = figsize_pixels)
    ax3 = CairoMakie.Axis(
        fig3[1, 1];
        xlabel = "x (m)", ylabel = "y (m)",
        title = "Vy (cm/yr x 1e4) Basic Grid",
        titlesize = 25,
        xlabelsize = 25, ylabelsize = 25,
        xticklabelsize = 20, yticklabelsize = 20,
        xticks = minimum(gridx_b):500:maximum(gridx_b),
        yticks = minimum(gridy_b):500:maximum(gridy_b),
    )
    CairoMakie.xlims!(ax3, minimum(gridx_b), maximum(gridx_b))
    CairoMakie.ylims!(ax3, minimum(gridy_b), maximum(gridy_b))
    hm3 = CairoMakie.heatmap!(
        ax3, gridx_b, gridy_b, permutedims(vyb_cm_yr) .* 1e4;
        colormap = :viridis,
    )
    CairoMakie.Colorbar(fig3[1, 2], hm3)
    plot_basic_grid_lines(ax3, model)
    plot_marker_overlay(ax3, model)
    CairoMakie.save("vy_basic_cm_yr_$(ntimestep_str).png", fig3)

    return nothing
end

function plot_viscoplastic_viscosity_etan0_etas0(model::ModelData, msg::String)::Nothing
    println("DEBUG: plotting viscoplastic viscosity etas0 and etan0 from ", msg)
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    ntimestep_str = string(ntimestep)
    gridx_b = model.grids.arrays.basic.gridx_b.array ./ 1000
    gridy_b = model.grids.arrays.basic.gridy_b.array ./ 1000
    gridx_pr = model.grids.arrays.pressure.gridx_pr.array ./ 1000
    gridy_pr = model.grids.arrays.pressure.gridy_pr.array ./ 1000
    etan0 = model.stokes_continuity.arrays.viscosity.etan0.array
    etas0 = model.stokes_continuity.arrays.viscosity.etas0.array

    ArrayUtils.print_min_max("etan0", etan0)
    ArrayUtils.print_min_max("etas0", etas0)
    
    dpi = 150
    figsize = (20, 5)
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)

    xspacing = 20
    yspacing = 20

    fig = CairoMakie.Figure(size = figsize_pixels)
    ax = CairoMakie.Axis(
        fig[1, 1];
        xlabel = "x (km)", ylabel = "y (km)",
        title = "Log10 Shear Viscosity (Pa.s)",
        titlesize = 25,
        xlabelsize = 25, ylabelsize = 25,
        xticklabelsize = 20, yticklabelsize = 20,
        xticks = minimum(gridx_b):xspacing:maximum(gridx_b),
        yticks = minimum(gridy_b):yspacing:maximum(gridy_b),
    )
    CairoMakie.xlims!(ax, minimum(gridx_b), maximum(gridx_b))
    CairoMakie.ylims!(ax, minimum(gridy_b), maximum(gridy_b))
    ax.yreversed = true
    hm = CairoMakie.heatmap!(
        ax, gridx_b, gridy_b, permutedims(log10.(etas0));
        colormap = :viridis,
    )
    CairoMakie.Colorbar(fig[1, 2], hm)
    CairoMakie.save("etas_2d_grid_$(ntimestep_str).png", fig)

    fig2 = CairoMakie.Figure(size = figsize_pixels)
    ax2 = CairoMakie.Axis(
        fig2[1, 1];
        xlabel = "x (km)", ylabel = "y (km)",
        title = "Log10 Normal Viscosity (Pa.s)",
        titlesize = 25,
        xlabelsize = 25, ylabelsize = 25,
        xticklabelsize = 20, yticklabelsize = 20,
        xticks = minimum(gridx_b):xspacing:maximum(gridx_b),
        yticks = minimum(gridy_b):yspacing:maximum(gridy_b),
    )
    CairoMakie.xlims!(ax2, minimum(gridx_b), maximum(gridx_b))
    CairoMakie.ylims!(ax2, minimum(gridy_b), maximum(gridy_b))
    ax2.yreversed = true
    hm2 = CairoMakie.heatmap!(
        ax2, gridx_pr, gridy_pr, permutedims(log10.(etan0));
        colormap = :viridis,
    )
    CairoMakie.Colorbar(fig2[1, 2], hm2)
    CairoMakie.save("etan_2d_grid_$(ntimestep_str).png", fig2)

    return nothing
end

const PLOT_MARKER_FRIC_MAX = 200_000  # max markers to plot; beyond this we downsample

"""
    plot_marker_fric(model::ModelData; filename::String="marker_fric_debug.png")

Scatter plot of `marker_fric` using **CairoMakie** (`scatter!` + **Colorbar**).
Positions in km; colormap `:turbo`; **colorrange** is fixed **0.0–0.5** (friction
coefficient). Decimates beyond `PLOT_MARKER_FRIC_MAX`.
"""
function plot_marker_fric(
    model::ModelData;
    filename::String="marker_fric_debug.png",
)::Nothing
    marker_x = model.markers.arrays.location.marker_x.array ./ 1000.0   # m -> km
    marker_y = model.markers.arrays.location.marker_y.array ./ 1000.0   # m -> km
    marker_fric = copy(model.markers.arrays.rheology.marker_fric.array)
    println("plot_marker_fric: marker_fric min: ", minimum(marker_fric), " max: ", maximum(marker_fric))
    colorrange_lo = 0.0
    colorrange_hi = 0.5
    cbar_ticks = collect(0.0:0.1:0.5)
    dpi = 150
    fig = CairoMakie.Figure(size=(round(Int, 10 * dpi), round(Int, 6 * dpi)))
    x_min, x_max = extrema(marker_x)
    xtick_positions = (floor(x_min / 25) * 25):25.0:(ceil(x_max / 25) * 25)
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel="x (km)",
        ylabel="y (km)",
        title="Marker friction coefficient (marker_fric)",
        aspect=CairoMakie.DataAspect(),
        yreversed=true,
        xticks=xtick_positions,
    )
    sc = CairoMakie.scatter!(
        ax,
        marker_x,
        marker_y;
        color=marker_fric,
        colormap=:turbo,
        colorrange=(colorrange_lo, colorrange_hi),
        markersize=2.0,
        marker=:circle,
        strokewidth=0.0,
        strokecolor=:black,
        rasterize=true,
    )
    cb = CairoMakie.Colorbar(fig[1, 2], sc, label="friction coef.")
    cb.ticks = cbar_ticks
    cb.labelsize = 14
    cb.ticklabelsize = 11
    CairoMakie.save(filename, fig)
    println(">> Saved marker_fric plot to $(filename)")
    return nothing
end

function plot_marker_overlay(ax::CairoMakie.Axis, model::ModelData)::Nothing
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    CairoMakie.scatter!(
        ax, marker_x, marker_y;
        color = :red, markersize = 1.0, strokewidth = 0.0,
    )
    return nothing
end

function plot_basic_grid_lines(ax::CairoMakie.Axis, model::ModelData)::Nothing
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array
    CairoMakie.vlines!(ax, gridx_b; color = :black, linewidth = 0.75)
    CairoMakie.hlines!(ax, gridy_b; color = :black, linewidth = 0.75)
    return nothing
end

end # module