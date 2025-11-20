module Strain

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_sediment_material_id
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_solidified_basalt_material_id
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_sticky_material_ids
import ...PlotMarkerArraysManager: PlotMarkerArrays
import ...Get: get_color_bar_axis_fractions
import ...Get: get_marker_coordinates
import ...MarkerColormapManager: MarkerColorMap
import ...FilterPlot: FilterPlotData
import ...FilterPlot: add_contour_description!
import ...FilterPlot: plot_filtered_marker_scalars_for_strain
import ...Ticks: get_colorbar_ticks
import ....PlotParametersManager: PlotParameters
import ....PlotParametersManager: update_plot_counter!
import ....PlotParametersManager.PlotViewManager: get_active_dimensions
import ....PlotDtypes: AxesType
import .....Scatter: plot_scatter

function plot_filtered_strain(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    axes::AxesType
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_plastic_strain
    plot_contours = parameters.marker_plot_params.plot_plastic_strain_contours

    marker_scalar_array = marker_arrays.marker_strain_plastic
    label = "Plastic Strain"

    strain_min = parameters.marker_plot_params.strain_min
    strain_max = parameters.marker_plot_params.strain_max
    min_and_max = (strain_min, strain_max)

    contour_interval = parameters.marker_plot_params.strain_contour_interval

    cmap_name = parameters.marker_plot_params.strain_cmap

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
        custom_cmap=nothing
    )

    sticky_matids = get_sticky_material_ids(materials)
    sed_matid = get_sediment_material_id(materials)
    basalt_matid = get_solidified_basalt_material_id(materials)

    add_contour_description!(filter_plot_data)

    if plot_marker_scalars == 1 || plot_contours == 1
        println(">> Plotting strain")
        plot_filtered_marker_scalars_for_strain(
            filter_plot_data, sticky_matids, sed_matid, basalt_matid)
    end

    return nothing
end

function plot_strain(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    colormap::MarkerColorMap,
    axes::AxesType
)::Nothing
    x_array_km, y_array_km = get_marker_coordinates(marker_arrays)
    color_array = copy(marker_arrays.marker_strain_plastic)
    marker_size = parameters.marker_plot_params.marker_size

    color_map = :None
    smax = 5.0
    nlevels = 5

    ticks = get_colorbar_ticks(smax, nlevels)

    n_bin = colormap.n_bin
    min_value = 0.0
    max_value = Float64(n_bin - 1)
    order_number = parameters.plot_counter
    
    plot_scatter(
        axes, color_map, x_array_km, y_array_km, color_array,
        marker_size, min_value, max_value, ticks;
        label="Plastic Strain",
        order_number=order_number,
        plot_dimensions=get_active_dimensions(parameters.view)
    )
    
    update_plot_counter!(parameters)
    
    return nothing
end

end # module
