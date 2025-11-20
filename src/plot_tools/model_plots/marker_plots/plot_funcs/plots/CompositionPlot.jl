module CompositionPlot

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: make_custom_labels
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: make_color_map_comp
import EarthBox.MaterialColorsContainer: get_colorbar_ticks_for_material_colors
import ...PlotMarkerArraysManager: PlotMarkerArrays
import ...MarkerColormapManager: MarkerColorMap
import ...Get: get_color_bar_axis_fractions
import ...Get: get_marker_coordinates
import ....PlotParametersManager: PlotParameters
import ....PlotParametersManager: update_plot_counter!
import ....PlotParametersManager: get_colorbar_order
import ....PlotParametersManager.PlotViewManager: get_active_dimensions
import ....PlotDtypes: AxesType
import .....Scatter: plot_scatter

function plot_composition(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    colormap::MarkerColorMap,
    axes::AxesType
)::Nothing
    println(">> Plotting composition")
    x_array_km, y_array_km = get_marker_coordinates(marker_arrays)
    color_array = copy(marker_arrays.marker_matid)

    marker_size = parameters.marker_plot_params.marker_size
    color_map = colormap.cm
    decimation_factor = parameters.marker_plot_params.decimation_factor

    n_bin = colormap.n_bin
    ticks = get_colorbar_ticks_for_material_colors(n_bin)

    min_value = 1.0 - 0.5
    max_value = Float64(n_bin) + 0.5

    order_number = get_colorbar_order(parameters)

    custom_labels = make_custom_labels(materials)

    colorbar_labels_fontsize = parameters.fonts.colorbar_labels_fontsize
    colorbar_ticks_fontsize = parameters.fonts.colorbar_ticks_fontsize
    
    plot_scatter(
        axes,
        color_map,
        x_array_km,
        y_array_km,
        color_array,
        marker_size,
        min_value,
        max_value,
        ticks;
        label="Material ID",
        order_number=order_number,
        custom_labels=custom_labels,
        colorbar_labels_fontsize=colorbar_labels_fontsize,
        colorbar_ticks_fontsize=colorbar_ticks_fontsize,
        decimation_factor=decimation_factor,
        plot_dimensions=get_active_dimensions(parameters.view),
        apply_domain_filter=true
    )

    update_plot_counter!(parameters)
    
    return nothing
end

end # module
