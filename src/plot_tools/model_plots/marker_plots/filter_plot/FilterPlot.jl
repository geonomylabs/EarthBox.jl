module FilterPlot

using Printf

import ..PlotMarkerArraysManager: PlotMarkerArrays
import ...PlotParametersManager: PlotParameters
import ...PlotParametersManager: PlotContoursManager
import ...PlotParametersManager: PlotViewManager
import ...PlotParametersManager: get_colorbar_order, update_plot_counter!
import ...PlotDtypes: AxesType
import ....Scatter: plot_scatter, make_grids_from_scatter
import ....Filters: filter_markers_based_on_minimum_value, 
                   filter_markers_for_strain_plot, 
                   filter_markers_for_age_plot
import ..Ticks: get_colorbar_ticks_general
import ..Get: get_color_bar_axis_fractions

mutable struct FilterPlotData
    parameters::Union{PlotParameters, Nothing}
    marker_arrays::Union{PlotMarkerArrays, Nothing}
    marker_scalar_array::Union{Vector{Float64}, Nothing}
    min_and_max::Union{Tuple{Float64, Float64}, Nothing}
    contour_interval::Union{Float64, Nothing}
    axes::Union{AxesType, Nothing}
    plot_marker_scalars::Union{Int, Nothing}
    plot_contours::Union{Int, Nothing}
    label::Union{String, Nothing}
    cmap_name::Union{String, Nothing}
    custom_cmap::Union{Any, Nothing}
    contour_type::Union{String, Nothing}
    discontinuous::Union{Bool, Nothing}
    matids_to_keep::Vector{Int16}
    number_format::String
    rightside_up::Bool

    function FilterPlotData(;
        parameters::Union{PlotParameters, Nothing} = nothing,
        marker_arrays::Union{PlotMarkerArrays, Nothing} = nothing,
        marker_scalar_array::Union{Vector{Float64}, Nothing} = nothing,
        min_and_max::Union{Tuple{Float64, Float64}, Nothing} = nothing,
        contour_interval::Union{Float64, Nothing} = nothing,
        axes::Union{AxesType, Nothing} = nothing,
        plot_marker_scalars::Union{Int, Nothing} = nothing,
        plot_contours::Union{Int, Nothing} = nothing,
        label::Union{String, Nothing} = nothing,
        cmap_name::Union{String, Nothing} = nothing,
        custom_cmap::Union{Any, Nothing} = nothing,
        contour_type::Union{String, Nothing} = nothing,
        discontinuous::Union{Bool, Nothing} = nothing,
        matids_to_keep::Vector{Int16} = Int16[],
        number_format::String = "%6.2f",
        rightside_up::Bool = true
    )
        return new(
            parameters, marker_arrays, marker_scalar_array, min_and_max,
            contour_interval, axes, plot_marker_scalars, plot_contours,
            label, cmap_name, custom_cmap, contour_type, discontinuous,
            matids_to_keep, number_format, rightside_up
        )
    end
end

function activate_discontinuous!(data::FilterPlotData)::Nothing
    data.discontinuous = true
    return nothing
end

function get_contour_color(data::FilterPlotData)::String
    if data.contour_type == "meltfrac"
        line_color = data.parameters.marker_plot_params.meltfrac_contour_color
    else
        line_color = data.parameters.marker_plot_params.contour_line_color
    end
    return line_color
end

function add_contour_description!(data::FilterPlotData)::Nothing
    if data.plot_contours == 1
        color = get_contour_color(data)
        label = data.label
        contour_interval = data.contour_interval
        description = "    " * color * " : " * label * " : CI=" * 
                     string(contour_interval)
        PlotContoursManager.add_to_contour_description!(
            data.parameters.contours, description)
    end
    return nothing
end

function plot_filtered_marker_scalars_based_on_minimum(
    filter_plot_data::FilterPlotData
)::Nothing
    scalar_min = filter_plot_data.min_and_max[1]
    matids_to_keep = filter_plot_data.matids_to_keep

    (
        marker_x_km_filtered,
        marker_y_km_filtered,
        marker_scalar_array_filtered,
        _marker_matid_filtered
    ) = filter_markers_based_on_minimum_value(
        scalar_min,
        filter_plot_data.marker_arrays.marker_x_km,
        filter_plot_data.marker_arrays.marker_y_km,
        filter_plot_data.marker_scalar_array,
        filter_plot_data.marker_arrays.marker_matid,
        matids_to_keep
    )
    if length(marker_scalar_array_filtered) == 0
        @printf("   >> No markers found for label = %s\n", filter_plot_data.label)
        # Create a dummy marker outside the plot domain:
        marker_x_km_filtered = [1.0e9]  # very large x, outside typical domain
        marker_y_km_filtered = [1.0e9]  # very large y, outside typical domain
        marker_scalar_array_filtered = [0.0]  # use a dummy scalar value
    end

    #if length(marker_scalar_array_filtered) > 0
    make_filter_plots(
        marker_x_km_filtered, marker_y_km_filtered,
        marker_scalar_array_filtered, filter_plot_data
    )
    #else
    #    @printf("   >> No markers found for label = %s\n", filter_plot_data.label)
    #end

    return nothing
end

function plot_filtered_marker_scalars_for_strain(
    filter_plot_data::FilterPlotData,
    sticky_matids::Tuple{Int16, Int16},
    sed_matid::Int16,
    basalt_matid::Int16
)::Nothing
    scalar_min = filter_plot_data.min_and_max[1]

    (
        marker_x_km_filtered,
        marker_y_km_filtered,
        marker_scalar_array_filtered,
        _marker_matid_filtered
    ) = filter_markers_for_strain_plot(
        scalar_min,
        filter_plot_data.marker_arrays.marker_x_km,
        filter_plot_data.marker_arrays.marker_y_km,
        filter_plot_data.marker_arrays.marker_pfailure,
        filter_plot_data.marker_scalar_array,
        filter_plot_data.marker_arrays.marker_matid,
        sticky_matids,
        sed_matid,
        basalt_matid
    )

    if length(marker_scalar_array_filtered) == 0
        @printf("   >> No markers found for label = %s\n", filter_plot_data.label)
        # Create a dummy marker outside the plot domain:
        marker_x_km_filtered = [1.0e9]  # very large x, outside typical domain
        marker_y_km_filtered = [1.0e9]  # very large y, outside typical domain
        marker_scalar_array_filtered = [0.0]  # use a dummy scalar value
    end

    #if length(marker_scalar_array_filtered) > 0
    make_filter_plots(
        marker_x_km_filtered, marker_y_km_filtered,
        marker_scalar_array_filtered, filter_plot_data
    )
    #else
    #    @printf("   >> No markers found for label = %s\n", filter_plot_data.label)
    #end

    return nothing
end

function plot_filtered_marker_scalars_for_age(
    filter_plot_data::FilterPlotData
)::Nothing
    (
        marker_x_km_filtered,
        marker_y_km_filtered,
        marker_scalar_array_filtered,
        _marker_matid_filtered
    ) = filter_markers_for_age_plot(
        filter_plot_data.marker_arrays.marker_x_km,
        filter_plot_data.marker_arrays.marker_y_km,
        filter_plot_data.marker_scalar_array,
        filter_plot_data.marker_arrays.marker_matid,
        filter_plot_data.matids_to_keep
    )

    if length(marker_scalar_array_filtered) == 0
        @printf("   >> No markers found for label = %s\n", filter_plot_data.label)
        # Create a dummy marker outside the plot domain:
        marker_x_km_filtered = [1.0e9]  # very large x, outside typical domain
        marker_y_km_filtered = [1.0e9]  # very large y, outside typical domain
        marker_scalar_array_filtered = [0.0]  # use a dummy scalar value
    end

    make_filter_plots(
        marker_x_km_filtered, marker_y_km_filtered,
        marker_scalar_array_filtered, filter_plot_data
    )

    return nothing
end

function make_filter_plots(
    marker_x_km_filtered::Vector{Float64},
    marker_y_km_filtered::Vector{Float64},
    marker_scalar_array_filtered::Vector{Float64},
    filter_plot_data::FilterPlotData
)::Nothing
    scalar_min = filter_plot_data.min_and_max[1]
    scalar_max = filter_plot_data.min_and_max[2]
    parameters = filter_plot_data.parameters

    marker_size = parameters.marker_plot_params.marker_size

    if filter_plot_data.plot_marker_scalars == 1
        color_bar_ticks = get_colorbar_ticks_general(
            scalar_min, scalar_max, filter_plot_data.contour_interval)

        order_number = get_colorbar_order(parameters)

        colorbar_labels_fontsize = parameters.fonts.colorbar_labels_fontsize
        colorbar_ticks_fontsize = parameters.fonts.colorbar_ticks_fontsize
        
        plot_scatter(
            filter_plot_data.axes,
            filter_plot_data.cmap_name,
            marker_x_km_filtered,
            marker_y_km_filtered,
            marker_scalar_array_filtered,
            marker_size,
            scalar_min,
            scalar_max,
            color_bar_ticks,
            label=filter_plot_data.label,
            custom_cmap=filter_plot_data.custom_cmap,
            order_number=order_number,
            colorbar_labels_fontsize=colorbar_labels_fontsize,
            colorbar_ticks_fontsize=colorbar_ticks_fontsize,
            plot_dimensions=PlotViewManager.get_active_dimensions(parameters.view),
            apply_domain_filter=true
        )
        update_plot_counter!(parameters)
    end

    if filter_plot_data.plot_contours == 1
        decimation_factor = parameters.marker_plot_params.decimation_factor_scatter_overlay
        nx = parameters.marker_plot_params.nx_contour_grid
        ny = parameters.marker_plot_params.ny_contour_grid

        grid_x, grid_y, grid_scalars = make_grids_from_scatter(
            marker_x_km_filtered, marker_y_km_filtered,
            marker_scalar_array_filtered; nx=nx, ny=ny,
            decimation_factor=decimation_factor
        )
        make_marker_scalar_contours(
            grid_x, grid_y, grid_scalars, filter_plot_data
        )
    end

    return nothing
end

function make_marker_scalar_contours(
    grid_x::Vector{Float64},
    grid_y::Vector{Float64},
    grid_scalars::Matrix{Float64},
    filter_plot_data::FilterPlotData
)::Nothing
    parameters = filter_plot_data.parameters
    PlotContoursManager.activate_contours!(parameters.contours)

    marker_plot_params = parameters.marker_plot_params

    parameters.contours.iplot_contour_labels = marker_plot_params.plot_contour_labels
    parameters.contours.contour_interval = filter_plot_data.contour_interval

    parameters.contours.excluded_vals = [0.0]

    min_and_max = filter_plot_data.min_and_max
    value_min = min_and_max[1]
    value_max = min_and_max[2]

    PlotContoursManager.update_contour_levels!(parameters.contours, value_min, value_max)
    PlotContoursManager.update_linewidths!(parameters.contours,
        value_min, value_max,
        marker_plot_params.contour_line_width
    )

    if filter_plot_data.contour_type == "meltfrac"
        line_color = marker_plot_params.meltfrac_contour_color
    else
        line_color = marker_plot_params.contour_line_color
    end

    contour_label_fontsize = parameters.fonts.contour_label_fontsize
    PlotContoursManager.plot_contours!(
        filter_plot_data.axes, 
        parameters.contours,
        grid_x, grid_y, grid_scalars,
        color=Symbol(line_color), 
        labelsize=contour_label_fontsize,
        number_format=filter_plot_data.number_format
    )

    return nothing
end

end # module
