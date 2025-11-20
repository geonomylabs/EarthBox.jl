module PlasticStrainRate

import ...PlotMarkerArraysManager: PlotMarkerArrays
import ...FilterPlot: FilterPlotData
import ...FilterPlot: add_contour_description!
import ...FilterPlot: plot_filtered_marker_scalars_based_on_minimum
import ....PlotParametersManager: PlotParameters
import ....PlotParametersManager.PlotViewManager: get_active_dimensions
import ....PlotDtypes: AxesType
import .....ColorMaps: make_transparent_single_color_colormap

function plot_filtered_plastic_strain_rate(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    axes::AxesType;
    use_log10::Bool = false
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_plastic_strain_rate
    plot_contours = parameters.marker_plot_params.plot_plastic_strain_rate_contours
    contour_interval = parameters.marker_plot_params.strain_rate_contour_interval
    strain_rate_min = parameters.marker_plot_params.strain_rate_min
    strain_rate_max = parameters.marker_plot_params.strain_rate_max

    if use_log10
        marker_scalar_array = marker_arrays.marker_strain_rate_plastic
        label = "log10(Plastic Strain Rate) (1/s)"
    else
        marker_scalar_array = 10.0 .^ marker_arrays.marker_strain_rate_plastic
        label = "Plastic Strain Rate (1/s)"
        strain_rate_min = 0.0
        strain_rate_max = 10.0^strain_rate_max
        contour_interval = (strain_rate_max - strain_rate_min) / 10.0
    end

    min_and_max = (strain_rate_min, strain_rate_max)

    cmap_name = parameters.marker_plot_params.strain_rate_cmap

    custom_cmap = make_transparent_single_color_colormap((1.0, 0.0, 0.0), max_alpha=1.0)

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
        println(">> Plotting plastic strain rate")
        plot_filtered_marker_scalars_based_on_minimum(filter_plot_data)
    end

    return nothing
end

end # module
