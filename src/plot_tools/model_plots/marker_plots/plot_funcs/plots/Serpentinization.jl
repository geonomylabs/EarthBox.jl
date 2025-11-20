module Serpentinization

import ...PlotMarkerArraysManager: PlotMarkerArrays
import ....PlotParametersManager: PlotParameters
import ....PlotDtypes: AxesType
import ...FilterPlot: FilterPlotData
import ...FilterPlot: add_contour_description!
import ...FilterPlot: plot_filtered_marker_scalars_based_on_minimum
import .....ColorMaps: make_transparent_single_color_colormap

""" Plot filtered serpentinization ratio for markers.

Plots serpentinization ratio markers with optional contours based on the 
provided parameters and marker arrays. Markers are filtered using a minimum 
serpentinization value.

Inputs
------
- parameters::PlotParameters: Plot configuration parameters
- marker_arrays::PlotMarkerArrays: Marker data arrays
- axes::AxesType: Plot axes object
"""
function plot_filtered_serpentinization(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    axes::AxesType
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_serpentinization
    plot_contours = parameters.marker_plot_params.plot_serpentinization_contours

    marker_scalar_array = marker_arrays.marker_serpentinization
    label = "Serpentinization Ratio"

    serpentinization_min = parameters.marker_plot_params.serpentinization_min
    serpentinization_max = parameters.marker_plot_params.serpentinization_max
    min_and_max = (serpentinization_min, serpentinization_max)

    contour_interval = parameters.marker_plot_params.serpentinization_contour_interval

    custom_cmap = make_transparent_single_color_colormap(
        (17.0/255.0, 107.0/255.0, 0.0))

    cmap_name = parameters.marker_plot_params.serpentinization_cmap

    matids_to_keep = Vector{Int16}(undef, 1)
    matids_to_keep[1] = -1
    
    filter_plot_data = FilterPlotData(
        parameters=parameters,
        marker_arrays=marker_arrays,
        marker_scalar_array=marker_scalar_array,
        min_and_max=min_and_max,
        contour_interval=contour_interval,
        axes=axes,
        plot_marker_scalars=plot_marker_scalars,
        plot_contours=plot_contours,
        label=label,
        cmap_name=cmap_name,
        custom_cmap=custom_cmap,
        matids_to_keep=matids_to_keep
    )

    add_contour_description!(filter_plot_data)

    if plot_marker_scalars == 1 || plot_contours == 1
        println(">> Plotting serpentinization")
        plot_filtered_marker_scalars_based_on_minimum(filter_plot_data)
    end

    return nothing
end

end # module
