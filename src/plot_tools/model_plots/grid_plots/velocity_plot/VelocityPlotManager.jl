module VelocityPlotManager

import CairoMakie
import CairoMakie: Colorbar
import EarthBox: JLDTools
import Statistics
import ..PlotVelocityArraysManager
import ...PlotTools
import ...PlotDtypes: PlotDictType, PlotParametersType
import ...PlotParametersManager
import ...PlotParametersManager: PlotConversionManager
import ...PlotParametersManager: PlotTimeManager

struct VelocityPlot
    dataname_x::String
    dataname_y::String
    dataname_mag::String
    plot_dict::PlotDictType
    parameter_group_name::String
    mainpath::String
    outpath::String
    parameters::PlotParametersManager.PlotParameters
    velocity_arrays::PlotVelocityArraysManager.PlotVelocityArrays
    figsize::Tuple{Float64, Float64}
end

function VelocityPlot(;
    dataname_x::String,
    dataname_y::String,
    dataname_mag::String,
    plot_dict::PlotDictType,
    parameter_group_name::String,
    mainpath::String,
    outpath::String
)::VelocityPlot
    return VelocityPlot(
        dataname_x, 
        dataname_y, 
        dataname_mag, 
        plot_dict, 
        parameter_group_name, 
        mainpath, 
        outpath, 
        PlotParametersManager.PlotParameters(plot_dict, mainpath, outpath), 
        PlotVelocityArraysManager.PlotVelocityArrays(), 
        (6.0, 4.0)
        )
end

function load_params(
    velocity_plot::VelocityPlot,
    plot_params::PlotParametersType
)::Nothing
    PlotParametersManager.update_scalar_plot_parameters!(
        velocity_plot.parameters, plot_params)
    return nothing
end

function make_2dgrid_plots_for_timestep_velocity(;
    velocity_plot::VelocityPlot,
    ioutput::Int64,
    decimation_factor::Int64=1,
    scale_factor::Float64=10.0
)::Nothing
    velocity_plot.parameters.time.ioutput = ioutput
    load_jld_velocity_information_and_set_arrays!(velocity_plot)
    println(">> Working on vector plot at time step $ioutput")
    make_2dvector_plot(
        velocity_plot=velocity_plot,
        decimation_factor=decimation_factor,
        scale_factor=scale_factor
    )
    return nothing
end

function make_2dvector_plot(;
    velocity_plot::VelocityPlot,
    decimation_factor::Int64=1,
    scale_factor::Float64=10.0
)::Nothing

    image_props = velocity_plot.parameters.image
    println(">> figsize: $(image_props.figsize)")

    gridx = velocity_plot.velocity_arrays.gridx
    gridy = velocity_plot.velocity_arrays.gridy
    Vectorx = velocity_plot.velocity_arrays.vectorx
    Vectory = velocity_plot.velocity_arrays.vectory
    
    gridx_decimated = gridx[1:decimation_factor:end]
    gridy_decimated = gridy[1:decimation_factor:end]
    vxb_cm_yr = Vectorx[1:decimation_factor:end, 1:decimation_factor:end]
    vyb_cm_yr = Vectory[1:decimation_factor:end, 1:decimation_factor:end]

    # Create coordinate arrays for arrow positions
    xnum = length(gridx_decimated)
    ynum = length(gridy_decimated)
    
    # Create meshgrid-like arrays for positions
    x_positions = zeros(Float64, ynum, xnum)
    y_positions = zeros(Float64, ynum, xnum)
    
    for j in 1:xnum
        for i in 1:ynum
            x_positions[i, j] = gridx_decimated[j]
            y_positions[i, j] = gridy_decimated[i]
        end
    end

    fig, axes_xy = PlotTools.InitializePlots.initialize_xy_plot(velocity_plot.parameters)

    velocity_magnitudes = calculate_velocity_magnitude(vxb_cm_yr, vyb_cm_yr)
    scaled_magnitudes = velocity_magnitudes .* scale_factor

    x_flat = vec(x_positions)
    y_flat = vec(y_positions)
    vx_flat = vec(vxb_cm_yr)
    vy_flat = vec(vyb_cm_yr)
    lengthscale_flat = vec(scaled_magnitudes)

    # Create the arrows plot with color based on magnitude
    arrows_plot = CairoMakie.arrows2d!(
        axes_xy, x_flat, y_flat, vx_flat, vy_flat,
        lengthscale=lengthscale_flat,
        color=vec(velocity_magnitudes),  # Color by velocity magnitude
        colormap=:rainbow,  # Choose a colormap
        shaftwidth=scale_factor/5.0,
        tipwidth=scale_factor,
        tiplength=scale_factor,
        alpha=0.8
    )

    # Add colorbar
    Colorbar(
        fig[1, 2], arrows_plot,
        label="Velocity Magnitude (cm/yr)",
        colorrange=(minimum(scaled_magnitudes), maximum(scaled_magnitudes)),
        )

    PlotTools.FinalizePlots.finalize_plot!(
        fig, axes_xy, velocity_plot.parameters, "Velocity_Field", 
        units="cm/yr", extension=".png"
        )
        
    return nothing
end

function plot_basic_grid_lines(
    ax::CairoMakie.Axis, gridx::Vector{Float64}, gridy::Vector{Float64}
)::Nothing
    # Plot vertical grid lines
    for x in gridx
        CairoMakie.lines!(ax, [x, x], [minimum(gridy), maximum(gridy)], 
            color=:black, linewidth=0.75)
    end
    # Plot horizontal grid lines 
    for y in gridy
        CairoMakie.lines!(ax, [minimum(gridx), maximum(gridx)], [y, y],
            color=:black, linewidth=0.75)
    end
    return nothing
end

function calculate_velocity_magnitude(
    vx::Matrix{Float64},
    vy::Matrix{Float64}
)::Matrix{Float64}
    ynum, xnum = size(vx)
    vmag = zeros(Float64, ynum, xnum)
    for i in 1:ynum
        for j in 1:xnum
            vyy = vy[i, j]
            vxx = vx[i, j]
            vmag[i, j] = sqrt(vxx*vxx + vyy*vyy)
        end
    end
    return vmag
end

function load_jld_velocity_information_and_set_arrays!(
    velocity_plot::VelocityPlot
)::Nothing
    jld_filename = get_jld_velocity_filename(velocity_plot)
    load_grid_and_time_using_velocity_x!(velocity_plot, jld_filename)
    veloc_x = get_velocity_component(velocity_plot, jld_filename, velocity_plot.dataname_x)
    veloc_y = get_velocity_component(velocity_plot, jld_filename, velocity_plot.dataname_y)
    veloc_mag = get_velocity_component(velocity_plot, jld_filename, velocity_plot.dataname_mag)
    PlotVelocityArraysManager.set_velocity_arrays!(
        velocity_plot.velocity_arrays, veloc_x, veloc_y, veloc_mag)
    return nothing
end

function load_grid_and_time_using_velocity_x!(
    vpt::VelocityPlot, 
    jld_filename::String
)::Nothing
    model_time, gridy, gridx, _, units_dict = JLDTools.get_jld_data(vpt.dataname_x, jld_filename)
    model_time, time_units = PlotConversionManager.convert_time_units(
        vpt.parameters.conversion, model_time, units_dict["time_units"])
    PlotTimeManager.set_plot_time_info!(vpt.parameters.time, model_time, time_units)
    gridx, gridy = PlotConversionManager.convert_grid_arrays_to_plot_units(
        vpt.parameters.conversion, gridx, gridy, units_dict["length_units"])
    PlotVelocityArraysManager.set_grid_arrays!(vpt.velocity_arrays, gridx, gridy)
    return nothing
end

function get_velocity_component(
    vpt::VelocityPlot,
    jld_filename::String,
    dataname::String
)::Matrix{Float64}
    _, _, _, ebarray2d, units_dict = JLDTools.get_jld_data(dataname, jld_filename)
    velocity_units = units_dict["scalar_units"]
    return PlotConversionManager.convert_velocity_array_units(
        vpt.parameters.conversion, velocity_units, ebarray2d)
end

function get_jld_velocity_filename(vpt::VelocityPlot)::String
    mainpath = vpt.parameters.paths.mainpath
    ioutput = vpt.parameters.time.ioutput
    return joinpath(mainpath, "vel_cmyr_" * JLDTools.intstr(ioutput) * ".jld")
end

end # module 