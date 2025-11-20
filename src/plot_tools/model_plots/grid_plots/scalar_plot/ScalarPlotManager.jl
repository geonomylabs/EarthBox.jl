module ScalarPlotManager

import CairoMakie
import EarthBox.PrintFuncs: print_warning
import EarthBox.ModelDataContainer: ModelData
import EarthBox.JLDTools: intstr, get_jld_data
import EarthBox.Arrays.ArrayUtils: check_limits
import EarthBox.GridFuncs: get_grid_info
import ...PlotDtypes: PlotDictType, PlotParametersType
import ...PlotParametersManager
import ...PlotParametersManager: PlotConversionManager
import ...PlotParametersManager: PlotTimeManager
import ...PlotParametersManager: PlotContoursManager
import ...PlotParametersManager: PlotColorBarManager
import ...PlotTools
import ..PlotScalarArraysManager
import ..PlotScalarArraysManager: PlotScalarArrays
import ..VelocityPlotManager
import ..PlotVelocityArraysManager: PlotVelocityArrays
import ..ScalarNamesManager: ScalarNames
import ..DataNames: ScalarH5Datanames, get_velocity_dataname_list
import ..ArrayLookupManager: ArrayLookup, get_array_info
import ..VelocityPlotManager: load_jld_velocity_information_and_set_arrays!

mutable struct ScalarPlot
    scalar_name::String
    dataname::String
    plot_dict::PlotDictType
    parameter_group_name::String
    mainpath::String
    outpath::String
    parameters::PlotParametersManager.PlotParameters
    scalar_arrays::PlotScalarArrays
    velocity_plot::VelocityPlotManager.VelocityPlot
    scalar_names::ScalarNames
end

function ScalarPlot(;
    scalar_name::String,
    dataname::String,
    plot_dict::PlotDictType,
    parameter_group_name::String,
    mainpath::String,
    outpath::String
)
    velocity_plot = VelocityPlotManager.VelocityPlot(
        dataname_x=ScalarH5Datanames().velx_cmyr,
        dataname_y=ScalarH5Datanames().vely_cmyr,
        dataname_mag=ScalarH5Datanames().velmag_cmyr,
        plot_dict=plot_dict,
        parameter_group_name=parameter_group_name,
        mainpath=mainpath,
        outpath=outpath
    )
    return ScalarPlot(
        scalar_name,
        dataname,
        plot_dict,
        parameter_group_name,
        mainpath,
        outpath,
        PlotParametersManager.PlotParameters(plot_dict, mainpath, outpath),
        PlotScalarArrays(),
        velocity_plot,
        ScalarNames()
    )
end

function set_ioutput!(scalar_plot::ScalarPlot, ioutput::Int64)::Nothing
    scalar_plot.parameters.time.ioutput = ioutput
    scalar_plot.velocity_plot.parameters.time.ioutput = ioutput
    return nothing
end

function load_params!(
    scalar_plot::ScalarPlot, 
    plot_parameters::PlotParametersType
)::Nothing
    PlotParametersManager.update_scalar_plot_parameters!(
        scalar_plot.parameters, plot_parameters)
end

function load_jld_information_and_set_arrays!(scalar_plot::ScalarPlot)
    scalar_datanames = [
        scalar_plot.velocity_plot.dataname_x,
        scalar_plot.velocity_plot.dataname_y,
        scalar_plot.velocity_plot.dataname_mag
    ]
    if !(scalar_plot.dataname in scalar_datanames)
        load_jld_field_information_and_set_arrays!(scalar_plot)
    else
        VelocityPlotManager.load_jld_velocity_information_and_set_arrays!(scalar_plot.velocity_plot)
        sync_scalar_plot_with_velocity_plot!(scalar_plot)
        PlotScalarArraysManager.set_grid_arrays!(
            scalar_plot.scalar_arrays,
            scalar_plot.velocity_plot.velocity_arrays.gridx,
            scalar_plot.velocity_plot.velocity_arrays.gridy
        )
        ebarray2d = if scalar_plot.dataname == scalar_plot.velocity_plot.dataname_x
            scalar_plot.velocity_plot.velocity_arrays.vectorx
        elseif scalar_plot.dataname == scalar_plot.velocity_plot.dataname_y
            scalar_plot.velocity_plot.velocity_arrays.vectory
        elseif scalar_plot.dataname == scalar_plot.velocity_plot.dataname_mag
            scalar_plot.velocity_plot.velocity_arrays.vectormag
        end
        
        PlotScalarArraysManager.set_scalar_arrays!(scalar_plot.scalar_arrays, ebarray2d)
    end
end

function sync_scalar_plot_with_velocity_plot!(scalar_plot::ScalarPlot)
    scalar_plot.parameters.time.plot_time = scalar_plot.velocity_plot.parameters.time.plot_time
    scalar_plot.parameters.time.plot_time_units = scalar_plot.velocity_plot.parameters.time.plot_time_units
    return nothing
end

function load_jld_field_information_and_set_arrays!(scalar_plot::ScalarPlot)
    h5_filename = get_h5_field_filename(scalar_plot)
    model_time, gridy, gridx, ebarray2d, units_dict = get_jld_data(scalar_plot.dataname, h5_filename)
    model_time, time_units = PlotConversionManager.convert_time_units(
        scalar_plot.parameters.conversion, model_time, units_dict["time_units"])
    PlotTimeManager.set_plot_time_info!(scalar_plot.parameters.time, model_time, time_units)

    gridx, gridy = PlotConversionManager.convert_grid_arrays_to_plot_units(
        scalar_plot.parameters.conversion, gridx, gridy, units_dict["length_units"])
    PlotScalarArraysManager.set_grid_arrays!(scalar_plot.scalar_arrays, gridx, gridy)

    units = units_dict["scalar_units"]
    # For viscosity and strainrate, the units are actually in
    # log10(Pa.s) and log10(1/s) if read from an H5 file
    if units == "Pa.s"
        units = "log10(Pa.s)"
    elseif units == "1/s"
        units = "log10(1/s)"
    end

    ebarray2d = convert_scalar_array_to_plot_units(scalar_plot, ebarray2d, units)
    PlotScalarArraysManager.set_scalar_arrays!(scalar_plot.scalar_arrays, ebarray2d)
end

function get_h5_field_filename(plot::ScalarPlot)
    mainpath = plot.parameters.paths.mainpath
    ioutput = plot.parameters.time.ioutput
    println("ioutput: ", ioutput)
    return joinpath(mainpath, "fields_" * intstr(ioutput) * ".jld")
end

function load_model_information_and_set_arrays!(
    scalar_plot::ScalarPlot, 
    model::ModelData
)::Nothing
    model_time = model.timestep.parameters.main_time_loop.timesum.value
    time_units = model.timestep.parameters.main_time_loop.timesum.units

    model_time, time_units = PlotConversionManager.convert_time_units(
        scalar_plot.parameters.conversion, model_time, time_units)
    PlotTimeManager.set_plot_time_info!(scalar_plot.parameters.time, model_time, time_units)

    array_lookup = ArrayLookup(model)
    ebarray2d, scalar_units, grid_type = get_array_info(array_lookup, scalar_plot.dataname)

    ebarray2d = convert_scalar_array_to_plot_units(scalar_plot, ebarray2d, scalar_units)
    units = scalar_plot.parameters.conversion.plot_units.active_units
    PlotScalarArraysManager.set_scalar_arrays!(scalar_plot.scalar_arrays, ebarray2d)
    
    if units === nothing
        units = "None"
    end
    check_limits(scalar_plot.scalar_name * " " * units, ebarray2d)

    gridx, gridy, length_units = get_grid_info(model.grids.arrays, grid_type)

    gridx, gridy = PlotConversionManager.convert_grid_arrays_to_plot_units(
        scalar_plot.parameters.conversion, gridx, gridy, length_units)
    PlotScalarArraysManager.set_grid_arrays!(scalar_plot.scalar_arrays, gridx, gridy)
    return nothing
end

function convert_scalar_array_to_plot_units(
    scalar_plot::ScalarPlot,
    ebarray2d::Matrix{Float64},
    scalar_units::String
)::Matrix{Float64}
    scalar_data_names = scalar_plot.scalar_arrays.scalar_data_names
    velocity_dataname_list = get_velocity_dataname_list(scalar_data_names)
    stress_names = [scalar_plot.scalar_names.shear_stress, scalar_plot.scalar_names.normal_stress]
    if scalar_plot.dataname in velocity_dataname_list
        println("scalar_units (velocity): ", scalar_units)
        ebarray2d = PlotConversionManager.convert_velocity_array_units(
            scalar_plot.parameters.conversion, scalar_units, ebarray2d)
    elseif scalar_plot.scalar_name == scalar_plot.scalar_names.temperature
        println("scalar_units (temperature): ", scalar_units)
        ebarray2d = PlotConversionManager.convert_temperature_array_units(
            scalar_plot.parameters.conversion, scalar_units, ebarray2d)
    elseif scalar_plot.scalar_name == scalar_plot.scalar_names.viscosity
        println("scalar_units (viscosity): ", scalar_units)
        ebarray2d = PlotConversionManager.convert_viscosity_array_units(
            scalar_plot.parameters.conversion, scalar_units, ebarray2d)
    elseif scalar_plot.scalar_name == scalar_plot.scalar_names.strainrate
        println("scalar_units (strainrate): ", scalar_units)
        ebarray2d = PlotConversionManager.convert_strainrate_array_units(
            scalar_plot.parameters.conversion, scalar_units, ebarray2d)
    elseif scalar_plot.scalar_name == scalar_plot.scalar_names.pressure
        println("scalar_units (pressure): ", scalar_units)
        ebarray2d = PlotConversionManager.convert_pressure_array_units(
            scalar_plot.parameters.conversion, scalar_units, ebarray2d)
    elseif scalar_plot.scalar_name in stress_names
        println("scalar_units (stress): ", scalar_units)
        ebarray2d = PlotConversionManager.convert_stress_array_units(
            scalar_plot.parameters.conversion, scalar_units, ebarray2d)
    end
    return ebarray2d
end

function make_2dgrid_plots_for_timestep(
    scalar_plot::ScalarPlot;
    ioutput::Union{Int64, Nothing}=nothing,
    model::Union{ModelData, Nothing}=nothing
)
    iplot = scalar_plot.parameters.options.iplot
    if iplot == 1
        if ioutput === nothing && model === nothing
            error("Both ioutput and model cannot be nothing.")
        end
        if ioutput !== nothing
            set_ioutput!(scalar_plot, ioutput)
            load_jld_information_and_set_arrays!(scalar_plot)
        elseif model !== nothing
            scalar_plot.parameters.time.ioutput = model.timestep.parameters.main_time_loop.ntimestep.value
            load_model_information_and_set_arrays!(scalar_plot, model)
        end
        
        PlotContoursManager.update_contour_levels!(
            scalar_plot.parameters.contours,
            scalar_plot.parameters.color_bar.minimum_value,
            scalar_plot.parameters.color_bar.maximum_value
        )
        PlotContoursManager.update_linewidths!(
            scalar_plot.parameters.contours,
            scalar_plot.parameters.color_bar.minimum_value,
            scalar_plot.parameters.color_bar.maximum_value
        )
        
        println(
            ">> Working on scalar plot $(scalar_plot.scalar_name) at time step " *
            "$(scalar_plot.parameters.time.ioutput)"
        )
        make_2dscalar_plot(scalar_plot)
    end
end

function make_2dscalar_plot(scalar_plot::ScalarPlot)
    grid_scalar = get_grid_scalar(scalar_plot)
    gridx = scalar_plot.scalar_arrays.gridx
    gridy = scalar_plot.scalar_arrays.gridy

    fig, axes_xy = PlotTools.InitializePlots.initialize_xy_plot(scalar_plot.parameters)
    color_map = scalar_plot.parameters.color_bar.color_map
    if color_map isa String
        color_map = Symbol(color_map)
    end
    # Note that with grid_scalar(i,j), i is in y-direction and j is in x-direction but
    # heatmap! expects gridx to be associated with the i-index and gridy to be associated
    # with the j-index. So a transpose is needed.
    #try
    #    eb_heatmap(axes_xy, gridx, gridy, grid_scalar, color_map, scalar_plot, true)
    #catch
        #print_warning("CairoMakie heatmap failed with interpolate=true probably because of irregular grid, trying interpolate=false", level=1)
    eb_heatmap(axes_xy, gridx, gridy, grid_scalar, color_map, scalar_plot, false)
    #end

    grid_plot_type = scalar_plot.parameters.options.grid_plot_type
    edgecolor = scalar_plot.parameters.options.edgecolor
    linewidth = scalar_plot.parameters.options.linewidth
    edgecolor = Symbol(edgecolor)
    if grid_plot_type == "mesh"
        for x in gridx
            CairoMakie.vlines!(axes_xy, x, color=(edgecolor, 0.3), linewidth=linewidth)
        end
        for y in gridy
            CairoMakie.hlines!(axes_xy, y, color=(edgecolor, 0.3), linewidth=linewidth)
        end
    end

    if scalar_plot.parameters.options.show_nodes == 1
        xcoors, ycoors = get_node_coordinates(scalar_plot)
        CairoMakie.scatter!(
            axes_xy, xcoors, ycoors, color=:black, markersize=1.0, 
            strokewidth=0.0, strokecolor=:transparent
            )
    end

    labelsize = scalar_plot.parameters.fonts.contour_label_fontsize
    PlotContoursManager.plot_contours!(
        axes_xy, scalar_plot.parameters.contours,
        gridx, gridy, grid_scalar, 
        labelsize=labelsize
    )
    active_units = scalar_plot.parameters.conversion.plot_units.active_units
    PlotColorBarManager.plot_colorbar!(
        fig, get_clims(scalar_plot), color_map, irow=1, icol=2, label=active_units)

    units = scalar_plot.parameters.conversion.plot_units.active_units
    extension = scalar_plot.parameters.image.extension
    PlotTools.FinalizePlots.finalize_plot!(
        fig, axes_xy, scalar_plot.parameters, scalar_plot.scalar_name, 
        units=units, extension=extension
        )
end

function eb_heatmap(
    axes_xy, gridx, gridy, grid_scalar, color_map, scalar_plot,
    interpolate::Bool = true
)
    CairoMakie.heatmap!(
        axes_xy, gridx, gridy, grid_scalar',
        colormap=color_map,
        colorrange=get_clims(scalar_plot),
        interpolate=interpolate # This does not work with irregular grids
        )
end

function get_grid_scalar(scalar_plot::ScalarPlot)::Matrix{Float64}
    return copy(scalar_plot.scalar_arrays.scalar)
end

function get_clims(scalar_plot::ScalarPlot)::Tuple{Float64, Float64}
    return (
        scalar_plot.parameters.color_bar.minimum_value, 
        scalar_plot.parameters.color_bar.maximum_value
    )
end

function get_node_coordinates(
    scalar_plot::ScalarPlot
)::Tuple{Vector{Float64}, Vector{Float64}}
    ynum = size(scalar_plot.scalar_arrays.gridy, 1)
    xnum = size(scalar_plot.scalar_arrays.gridx, 1)
    ncoors = xnum * xnum
    xcoors = zeros(Float64, ncoors)
    ycoors = zeros(Float64, ncoors)
    icount = 1
    
    for i in 1:xnum
        for j in 1:ynum
            xcoors[icount] = scalar_plot.scalar_arrays.gridx[i]
            ycoors[icount] = scalar_plot.scalar_arrays.gridy[j]
            icount += 1
        end
    end
    return xcoors, ycoors
end

function print_array_shapes(
    scalar_plot::ScalarPlot, 
    grid_scalar::Array{Float64}
)::Nothing
    println("y-grid shape: $(size(scalar_plot.scalar_arrays.gridy))")
    println("x-grid shape: $(size(scalar_plot.scalar_arrays.gridx))")
    println("scalar shape: $(size(grid_scalar))")
    return nothing
end

function print_array_dimensions(scalar_plot::ScalarPlot)::Nothing
    gridy = scalar_plot.scalar_arrays.gridy
    gridx = scalar_plot.scalar_arrays.gridx
    println("y-grid: min $(minimum(gridy)) max $(maximum(gridy))")
    println("x-grid: min $(minimum(gridx)) max $(maximum(gridx))")
    return nothing
end

end # module 