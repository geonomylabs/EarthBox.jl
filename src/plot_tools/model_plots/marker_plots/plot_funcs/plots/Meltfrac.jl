module Meltfrac

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_mantle_melting_matids
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_mantle_and_gabbroic_melting_matids
import ...PlotMarkerArraysManager: PlotMarkerArrays
import ....PlotParametersManager: PlotParameters
import ....PlotDtypes: AxesType
import ...FilterPlot: FilterPlotData
import ...FilterPlot: add_contour_description!
import ...FilterPlot: plot_filtered_marker_scalars_based_on_minimum

""" Plot marker meltfrac.

Plots melt fraction markers with optional contours based on the provided 
parameters and marker arrays. Filters markers to show either mantle melting 
or mantle and gabbroic melting material IDs.

Inputs
------
- parameters::PlotParameters: Plot configuration parameters
- marker_arrays::PlotMarkerArrays: Marker data arrays
- materials::Materials: Material properties container
- axes::AxesType: Plot axes object

Keyword Arguments
----------------
- use_gabbro_melting::Bool: Whether to include gabbroic melting materials 
  (default: false)
"""
function plot_filtered_meltfrac(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    axes::AxesType;
    use_gabbro_melting::Bool = false
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_meltfrac
    plot_contours = parameters.marker_plot_params.plot_meltfrac_contours
    number_format = parameters.marker_plot_params.meltfrac_number_format
    rightside_up = parameters.marker_plot_params.meltfrac_label_rightside_up
    marker_scalar_array = marker_arrays.marker_meltfrac
    label = "Melt Fraction"

    meltfrac_min = parameters.marker_plot_params.melt_fraction_min
    meltfrac_max = parameters.marker_plot_params.melt_fraction_max
    min_and_max = (meltfrac_min, meltfrac_max)

    contour_interval = parameters.marker_plot_params.melt_fraction_contour_interval

    cmap_name = parameters.marker_plot_params.meltfrac_cmap

    if use_gabbro_melting
        matids_to_keep = get_mantle_and_gabbroic_melting_matids(materials)
    else
        matids_to_keep = get_mantle_melting_matids(materials)
    end

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
        contour_type="meltfrac",
        matids_to_keep=matids_to_keep,
        number_format=number_format,
        rightside_up=rightside_up
    )

    add_contour_description!(filter_plot_data)

    if plot_marker_scalars == 1 || plot_contours == 1
        println(">> Plotting meltfrac")
        plot_filtered_marker_scalars_based_on_minimum(filter_plot_data)
    end

    return nothing
end

end # module
