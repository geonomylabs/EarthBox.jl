module DebugPlots

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Arrays: ArrayUtils
import EarthBox.ConversionFuncs: get_factor_cm_yr_to_m_s
import Plots

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

    p = Plots.heatmap(
        gridx_b, gridy_b, tk1, size=figsize_pixels,
        xlabel="x (km)", ylabel="y (km)", 
        title="Temperature (K)",
        color=:viridis, aspect_ratio=:auto,
        xlims=(minimum(gridx_b), maximum(gridx_b)),
        ylims=(minimum(gridy_b), maximum(gridy_b)),
        xticks=minimum(gridx_b):xspacing:maximum(gridx_b),
        yticks=minimum(gridy_b):yspacing:maximum(gridy_b),
        legendfontsize=25, guidefontsize=25,
        tickfontsize=20, titlefontsize=25,
    )
    Plots.plot!(p, yflip=true)
    #plot_basic_grid_lines(p, model)
    #plot_marker_overlay(p, model)
    Plots.savefig("tk1_2d_grid_$(ntimestep_str).png")
    
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
    
    p = Plots.quiver(
        x, y, quiver=(vx, vy), size=figsize_pixels,
        xlabel="x (m)", ylabel="y (m)",
        title="Basic Grid Velocity Field (Normalized)", aspect_ratio=:auto,
        xlims=(minimum(gridx_b), maximum(gridx_b)),
        ylims=(minimum(gridy_b), maximum(gridy_b)),
        xticks=minimum(gridx_b):500:maximum(gridx_b),
        yticks=minimum(gridy_b):500:maximum(gridy_b),
        legendfontsize=12, guidefontsize=12, tickfontsize=10, titlefontsize=15
    )
    Plots.savefig("basic_grid_velocity_field_$(ntimestep_str).png")

    p = Plots.heatmap(
        gridx_b, gridy_b, vxb_cm_yr.*1e4, size=figsize_pixels,
        dip=300, xlabel="x (m)", ylabel="y (m)", 
        title="Vx (cm/yr x 1e4) Basic Grid",
        color=:viridis, aspect_ratio=:auto,
        xlims=(minimum(gridx_b), maximum(gridx_b)),
        ylims=(minimum(gridy_b), maximum(gridy_b)),
        xticks=minimum(gridx_b):500:maximum(gridx_b),
        yticks=minimum(gridy_b):500:maximum(gridy_b),
        legendfontsize=25, guidefontsize=25,
        tickfontsize=20, titlefontsize=25,
    )
    plot_basic_grid_lines(p, model)
    plot_marker_overlay(p, model)
    Plots.savefig("vx_basic_cm_yr_$(ntimestep_str).png")

    p = Plots.heatmap(
        gridx_b, gridy_b, vyb_cm_yr.*1e4, size=figsize_pixels,
        dip=300, xlabel="x (m)", ylabel="y (m)", 
        title="Vy (cm/yr x 1e4) Basic Grid",
        color=:viridis, aspect_ratio=:auto,
        xlims=(minimum(gridx_b), maximum(gridx_b)),
        ylims=(minimum(gridy_b), maximum(gridy_b)),
        xticks=minimum(gridx_b):500:maximum(gridx_b),
        yticks=minimum(gridy_b):500:maximum(gridy_b),
        legendfontsize=25, guidefontsize=25,
        tickfontsize=20, titlefontsize=25,
    )
    plot_basic_grid_lines(p, model)
    plot_marker_overlay(p, model)
    Plots.savefig("vy_basic_cm_yr_$(ntimestep_str).png")

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

    p = Plots.heatmap(
        gridx_b, gridy_b, log10.(etas0), size=figsize_pixels,
        xlabel="x (km)", ylabel="y (km)", 
        title="Log10 Shear Viscosity (Pa.s)",
        color=:viridis, aspect_ratio=:auto,
        xlims=(minimum(gridx_b), maximum(gridx_b)),
        ylims=(minimum(gridy_b), maximum(gridy_b)),
        xticks=minimum(gridx_b):xspacing:maximum(gridx_b),
        yticks=minimum(gridy_b):yspacing:maximum(gridy_b),
        legendfontsize=25, guidefontsize=25,
        tickfontsize=20, titlefontsize=25,
    )
    Plots.plot!(p, yflip=true)
    #plot_basic_grid_lines(p, model)
    #plot_marker_overlay(p, model)
    Plots.savefig("etas_2d_grid_$(ntimestep_str).png")
    
    p = Plots.heatmap(
        gridx_pr, gridy_pr, log10.(etan0), size=figsize_pixels,
        xlabel="x (km)", ylabel="y (km)", 
        title="Log10 Normal Viscosity (Pa.s)",
        color=:viridis, aspect_ratio=:auto,
        xlims=(minimum(gridx_b), maximum(gridx_b)),
        ylims=(minimum(gridy_b), maximum(gridy_b)),
        xticks=minimum(gridx_b):xspacing:maximum(gridx_b),
        yticks=minimum(gridy_b):yspacing:maximum(gridy_b),
        legendfontsize=25, guidefontsize=25,
        tickfontsize=20, titlefontsize=25,
    )
    Plots.plot!(p, yflip=true)
    #plot_basic_grid_lines(p, model)
    #plot_marker_overlay(p, model)
    Plots.savefig("etan_2d_grid_$(ntimestep_str).png")

    return nothing
end

function plot_marker_overlay(p::Plots.Plot, model::ModelData)::Nothing
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    Plots.scatter!(
        p, marker_x, marker_y, 
        color=:red,
        markersize=0.25,
        label=nothing
    )
    return nothing
end

function plot_basic_grid_lines(p::Plots.Plot, model::ModelData)::Nothing
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array
    # Plot vertical grid lines
    for x in gridx_b
        Plots.plot!(p, [x, x], [minimum(gridy_b), maximum(gridy_b)], 
            color=:black, label=nothing, linewidth=0.75)
    end
    # Plot horizontal grid lines 
    for y in gridy_b
        Plots.plot!(p, [minimum(gridx_b), maximum(gridx_b)], [y, y],
            color=:black, label=nothing, linewidth=0.75)
    end
end

end # module