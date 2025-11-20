module PlasticFailure

import ...PlotMarkerArraysManager: PlotMarkerArrays
import ...MarkerColormapManager: MarkerColorMap
import ...Get: get_marker_coordinates
import ...Ticks: get_colorbar_ticks_for_composition_plot
import ....PlotParametersManager: PlotParameters
import ....PlotParametersManager: update_plot_counter!
import ....PlotParametersManager: get_colorbar_order
import ....PlotParametersManager.PlotViewManager: get_active_dimensions
import ....PlotDtypes: AxesType
import .....Scatter: plot_scatter

function plot_plastic(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    colormap::MarkerColorMap,
    axes::AxesType
)::Nothing
    println(">> Plotting plastic failure")
    x_array_km, y_array_km = get_marker_coordinates(marker_arrays)
    color_array = copy(marker_arrays.marker_pfailure)
    
    marker_size = parameters.marker_plot_params.marker_size
    # TODO: update this color map definition for CairoMakie. It was defined as
    # 'None' in the original code.
    color_map = "inferno"
    ticks = get_colorbar_ticks_for_composition_plot(colormap.n_bin)
    
    min_value = 0.0
    max_value = 1.0

    order_number = get_colorbar_order(parameters)
    
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
        order_number=order_number,
        colorbar_labels_fontsize=colorbar_labels_fontsize,
        colorbar_ticks_fontsize=colorbar_ticks_fontsize,
        plot_dimensions=get_active_dimensions(parameters.view)
    )

    update_plot_counter!(parameters)
    
    return nothing
end

end # module
