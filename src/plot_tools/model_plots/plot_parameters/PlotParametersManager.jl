module PlotParametersManager

include("parameters/PlotFontsManager.jl")
include("parameters/PlotUnitsManager.jl")
include("parameters/PlotColorBarManager.jl")
include("parameters/PlotViewManager.jl")
include("parameters/PlotTicksManager.jl")
include("parameters/PlotConversionManager.jl")
include("parameters/PlotContoursManager.jl")
include("parameters/PlotPathsManager.jl")
include("parameters/PlotTimeManager.jl")
include("parameters/PlotOptionsManager.jl")
include("parameters/ImagePropsManager.jl")
include("parameters/MarkerPlotParametersManager.jl")


import EarthBox.UnitConversion: UnitConversionData
import ..PlotDtypes: PlotDictType, PlotParametersType, AbstractPlotParameterGroup
import .PlotFontsManager: PlotFonts
import .PlotUnitsManager: PlotUnits
import .PlotColorBarManager: PlotColorBar
import .PlotConversionManager: PlotConversions
import .PlotContoursManager
import .PlotContoursManager: PlotContours
import .PlotPathsManager: PlotPaths
import .PlotTimeManager: PlotTime
import .PlotOptionsManager: PlotOptions
import .ImagePropsManager: ImageProps
import .PlotViewManager: PlotView, update_active_zoom!
import .PlotTicksManager: PlotTicks, update_ticks!
import .MarkerPlotParametersManager: MarkerPlotParameters

mutable struct PlotParameters
    plot_dict::PlotDictType
    options::PlotOptions
    time::PlotTime
    paths::PlotPaths
    contours::PlotContours
    color_bar::PlotColorBar
    conversion::PlotConversions
    image::ImageProps
    view::PlotView
    ticks::PlotTicks
    units::PlotUnits
    marker_plot_params::MarkerPlotParameters
    unit_conversion::UnitConversionData
    fonts::PlotFonts
    plot_counter::Int
    y_sealevel::Float64
    base_level_shift::Float64
    total_plots::Int
end

function PlotParameters(
    plot_dict::PlotDictType,
    mainpath::String,
    outpath::String
)::PlotParameters
    view = PlotView(plot_dict)
    update_active_zoom!(view)
    ticks = PlotTicks(plot_dict)
    update_ticks!(ticks, view)
    return PlotParameters(
        plot_dict, 
        PlotOptions(plot_dict), 
        PlotTime(), 
        PlotPaths(mainpath=mainpath, outpath=outpath),
        PlotContours(plot_dict), 
        PlotColorBar(plot_dict), 
        PlotConversions(plot_dict), 
        ImageProps(plot_dict), 
        view, 
        ticks, 
        PlotUnits(plot_dict), 
        MarkerPlotParameters(),
        UnitConversionData(), 
        PlotFonts(plot_dict), 
        0, 
        0.0, 
        0.0,
        0
    )
end

function set_y_sealevel!(params::PlotParameters, y_sealevel::Float64)::Nothing
    params.y_sealevel = y_sealevel
    return nothing
end

function set_base_level_shift!(params::PlotParameters, base_level_shift::Float64)::Nothing
    params.base_level_shift = base_level_shift
    return nothing
end

function get_colorbar_order(params::PlotParameters)::Int
    total_colorbars = count_number_of_colorbars_to_plot(params)
    plot_counter = params.plot_counter
    #return total_colorbars - 1 - plot_counter
    return 2 + plot_counter
end

function count_number_of_colorbars_to_plot(params::PlotParameters)::Int
    return (
        params.marker_plot_params.plot_plastic_strain
        + params.marker_plot_params.plot_plastic_strain_rate
        + params.marker_plot_params.plot_sediment_age
        + params.marker_plot_params.plot_volcanics_age
        + params.marker_plot_params.plot_intrusive_age
        + params.marker_plot_params.plot_meltfrac
        + params.marker_plot_params.plot_extracted_meltfrac
        + params.marker_plot_params.plot_extractable_meltfrac
        + params.marker_plot_params.plot_serpentinization
        + 1
        )
end

function update_plot_counter!(params::PlotParameters)::Nothing
    params.plot_counter += 1
    return nothing
end

function reset_plot_counter!(params::PlotParameters)::Nothing
    params.plot_counter = 0
    return nothing
end

function get_plot_units(
    params::PlotParameters
)::Tuple{Union{String, Nothing}, Union{String, Nothing}}
    return (
        params.conversion.plot_units.length_units,
        params.conversion.plot_units.time_units
    )
end

function set_plot_time!(params::PlotParameters, model_time::Float64)::Nothing
    params.time.tMyr = model_time
    return nothing
end

function update_scalar_plot_parameters!(
    params::PlotParameters,
    parameters::PlotParametersType
)::Nothing
    set_parameter_group_attributes!(params.contours, parameters)
    PlotContoursManager.update_excluded_values_list!(params.contours)
    set_parameter_group_attributes!(params.color_bar, parameters)
    set_parameter_group_attributes!(params.options, parameters)
    return nothing
end

function update_marker_plot_parameters!(params::PlotParameters)::Nothing
    set_parameter_group_attributes!(params.marker_plot_params, params.plot_dict["marker_plot"])
    return nothing
end

function set_parameter_group_attributes!(
    parameter_group::AbstractPlotParameterGroup,
    parameters::PlotParametersType
)::Nothing
    for field in fieldnames(typeof(parameter_group))
        #println("Working on field $field of parameter group")
        if String(field) in keys(parameters)
            # Check that the type matches before setting
            field_type = typeof(getfield(parameter_group, field))
            param_value = parameters[String(field)]
            if param_value isa Tuple
                param_value = collect(param_value)
            end
            if !isa(param_value, field_type)
                # println(">> Type mismatch for field $field with value $param_value: expected $field_type but got $(typeof(param_value)). Trying to convert...")
                try
                    param_value = convert(field_type, param_value)
                catch
                    error(">> Failed to convert value $param_value to type $field_type for field $field")
                end
            end
            #println("Setting field $field to $param_value")
            setfield!(parameter_group, field, param_value)
        end
    end
    return nothing
end

end # module 