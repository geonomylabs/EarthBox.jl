module MarkerAge

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_sediment_material_id
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_solidified_basalt_material_id
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_solidified_layered_gabbro_material_id
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_solidified_gabbro_material_id
import ...PlotMarkerArraysManager: PlotMarkerArrays
import ...Ticks: get_colorbar_ticks_general
import ....PlotParametersManager: PlotParameters
import ....PlotDtypes: AxesType
import .....ColorMaps: make_alternating_colormap, make_discontinuous_colormap_from_input_cmap
import ...FilterPlot: FilterPlotData
import ...FilterPlot: activate_discontinuous!
import ...FilterPlot: add_contour_description!
import ...FilterPlot: plot_filtered_marker_scalars_for_age

function plot_filtered_sediment_age(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    axes::AxesType
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_sediment_age
    plot_contours = parameters.marker_plot_params.plot_sediment_age_contours

    marker_scalar_array = marker_arrays.marker_age
    label = "Deposition Time (Myr)"

    age_min = parameters.marker_plot_params.age_min
    age_max = parameters.marker_plot_params.age_max
    min_and_max = (age_min, age_max)

    contour_interval = parameters.marker_plot_params.age_contour_interval

    matids_to_keep = Vector{Int16}(undef, 1)
    matids_to_keep[1] = get_sediment_material_id(materials)
    cmap_name = parameters.marker_plot_params.age_cmap

    custom_cmap = make_discontinuous_colormap_from_input_cmap(
        age_min, age_max, contour_interval, cmap_name)
    cmap_name = "None"

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
        println(">> Plotting sediment age")
        plot_filtered_marker_scalars_for_age(filter_plot_data)
    end

    return nothing
end

function plot_filtered_volcanics_age(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    axes::AxesType
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_volcanics_age
    plot_contours = parameters.marker_plot_params.plot_volcanics_age_contours
    use_alternating_colormap_volcanics = 
        parameters.marker_plot_params.use_alternating_colormap_volcanics
    marker_scalar_array = marker_arrays.marker_age
    label = "Flow Solidification Time (Myr)"

    age_min = parameters.marker_plot_params.age_min_volcanics
    age_max = parameters.marker_plot_params.age_max_volcanics
    min_and_max = (age_min, age_max)

    contour_interval = parameters.marker_plot_params.age_contour_interval_volcanics

    matids_to_keep = Vector{Int16}(undef, 1)
    matids_to_keep[1] = get_solidified_basalt_material_id(materials)

    if use_alternating_colormap_volcanics
        custom_cmap = make_alternating_colormap(
            age_min, age_max, contour_interval,
            (1.0, 0.855, 0.725), (1.0, 0.498, 0.314)
        )
        cmap_name = "None"
    else
        custom_cmap = make_discontinuous_colormap_from_input_cmap(
            age_min, age_max, contour_interval, cmap_name)
        cmap_name = "None"
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
        custom_cmap=custom_cmap,
        matids_to_keep=matids_to_keep
    )

    add_contour_description!(filter_plot_data)

    if plot_marker_scalars == 1 || plot_contours == 1
        println(">> Plotting volcanics age")
        plot_filtered_marker_scalars_for_age(filter_plot_data)
    end

    return nothing
end

function plot_filtered_intrusive_age(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    axes::AxesType
)::Nothing
    plot_marker_scalars = parameters.marker_plot_params.plot_intrusive_age
    plot_contours = parameters.marker_plot_params.plot_intrusive_age_contours
    use_alternating_colormap_intrusive = 
        parameters.marker_plot_params.use_alternating_colormap_intrusive

    marker_scalar_array = marker_arrays.marker_age
    label = "Intrusive Solidification Time (Myr)"

    age_min = parameters.marker_plot_params.age_min_intrusive
    age_max = parameters.marker_plot_params.age_max_intrusive
    min_and_max = (age_min, age_max)

    contour_interval = parameters.marker_plot_params.age_contour_interval_intrusive

    matids_to_keep = Vector{Int16}(undef, 2)
    matids_to_keep[1] = get_solidified_layered_gabbro_material_id(materials)
    matids_to_keep[2] = get_solidified_gabbro_material_id(materials)

    if use_alternating_colormap_intrusive
        custom_cmap = make_alternating_colormap(
            age_min, age_max, contour_interval,
            (0.8, 0.8, 0.8), (0.4, 0.4, 0.4)
        )
        cmap_name = "None"
    else
        custom_cmap = make_discontinuous_colormap_from_input_cmap(
            age_min, age_max, contour_interval, cmap_name)
        cmap_name = "None"
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
        custom_cmap=custom_cmap,
        matids_to_keep=matids_to_keep
    )

    add_contour_description!(filter_plot_data)

    if plot_marker_scalars == 1 || plot_contours == 1
        println(">> Plotting intrusive age")
        plot_filtered_marker_scalars_for_age(filter_plot_data)
    end

    return nothing
end

end # module
