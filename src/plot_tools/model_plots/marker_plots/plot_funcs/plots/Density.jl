module Density

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_mantle_melting_matids
import ...PlotMarkerArraysManager: PlotMarkerArrays
import ....PlotParametersManager: PlotParameters
import ....PlotDtypes: AxesType
import ...FilterPlot: FilterPlotData
import ...FilterPlot: add_contour_description!
import ...FilterPlot: plot_filtered_marker_scalars_based_on_minimum

""" Plot filtered density.

Plots density markers with optional contours based on the provided parameters 
and marker arrays.

Inputs
------
- parameters::PlotParameters: Plot configuration parameters
- marker_arrays::PlotMarkerArrays: Marker data arrays
- axes::AxesType: Plot axes object
"""
function plot_filtered_density(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    axes::AxesType
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_density
    plot_contours = parameters.marker_plot_params.plot_density_contours

    marker_scalar_array = marker_arrays.marker_rho
    label = "Density"

    density_min = parameters.marker_plot_params.density_min
    density_max = parameters.marker_plot_params.density_max
    min_and_max = (density_min, density_max)

    contour_interval = parameters.marker_plot_params.density_contour_interval
    cmap_name = parameters.marker_plot_params.density_cmap

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
        custom_cmap=nothing,
        matids_to_keep=matids_to_keep
    )

    add_contour_description!(filter_plot_data)

    if plot_marker_scalars == 1 || plot_contours == 1
        println(">> Plotting density")
        plot_filtered_marker_scalars_based_on_minimum(filter_plot_data)
    end

    return nothing
end

""" Plot filtered mantle density.

Plots mantle density markers with optional contours based on the provided 
parameters and marker arrays. Filters markers to only show those with 
mantle melting material IDs.

Inputs
------
- parameters::PlotParameters: Plot configuration parameters
- marker_arrays::PlotMarkerArrays: Marker data arrays
- materials::Materials: Material properties container
- axes::AxesType: Plot axes object
"""
function plot_filtered_mantle_density(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    axes::AxesType
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_mantle_density
    plot_contours = parameters.marker_plot_params.plot_mantle_density_contours

    marker_scalar_array = marker_arrays.marker_rho
    label = "Mantle Density (kg/m^3)"

    density_min = parameters.marker_plot_params.density_min
    density_max = parameters.marker_plot_params.density_max
    min_and_max = (density_min, density_max)

    contour_interval = parameters.marker_plot_params.density_contour_interval
    cmap_name = parameters.marker_plot_params.density_cmap

    matids_to_keep = get_mantle_melting_matids(materials)
    
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
        custom_cmap=nothing,
        matids_to_keep=matids_to_keep
    )

    add_contour_description!(filter_plot_data)

    if plot_marker_scalars == 1 || plot_contours == 1
        plot_filtered_marker_scalars_based_on_minimum(filter_plot_data)
    end

    return nothing
end

end # module
