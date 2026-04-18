module CompositionPlot

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.MaterialsContainer: make_custom_labels
import EarthBox.MaterialColorsContainer: get_colorbar_ticks_for_material_colors
import ...PlotMarkerArraysManager: PlotMarkerArrays
import ...MarkerColormapManager: MarkerColorMap
import ...Get: get_marker_coordinates
import ....PlotParametersManager: PlotParameters
import ....PlotParametersManager: update_plot_counter!
import ....PlotParametersManager: get_colorbar_order
import ....PlotParametersManager.PlotViewManager: get_active_dimensions
import ....PlotDtypes: AxesType
import .....Scatter: plot_scatter

const REDUNDANT_MAT_TYPE_PHRASES = [
    "Sticky Air",
    "Sticky Water",
    "Sediment",
    "Felsic Continental Crust",
    "Mafic Continental Crust",
    "Ultramafic Mantle",
    "Layered Gabbro",
]

function filter_redundant_type_from_labels(labels::Vector{String})::Vector{String}
    return map(labels) do label
        m = match(r"^(.*) \((.+)\)$", label)
        if isnothing(m)
            return label
        end
        mat_part = String(m.captures[1])
        type_part = String(m.captures[2])
        for phrase in REDUNDANT_MAT_TYPE_PHRASES
            type_part = replace(type_part, phrase => "")
        end
        type_part = strip(replace(type_part, r" +" => " "))
        isempty(type_part) ? mat_part : "$mat_part ($type_part)"
    end
end

function plot_composition(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    colormap::MarkerColorMap,
    axes::AxesType
)::Nothing
    println(">> Plotting composition")
    x_array_km, y_array_km = get_marker_coordinates(marker_arrays)
    n_bin = colormap.n_bin
    # Reversed categorical colormap: map matid so each material keeps its RGB after reverse!.
    color_array = Float64(n_bin) .+ 1.0 .- copy(marker_arrays.marker_matid)

    marker_size = parameters.marker_plot_params.marker_size
    color_map = colormap.cm
    decimation_factor = parameters.marker_plot_params.decimation_factor

    ticks = get_colorbar_ticks_for_material_colors(n_bin)

    min_value = 1.0 - 0.5
    max_value = Float64(n_bin) + 0.5

    order_number = get_colorbar_order(parameters)

    custom_labels = make_custom_labels(materials)

    custom_labels = filter_redundant_type_from_labels(custom_labels)

    if !parameters.marker_plot_params.show_composition_matid_labels
        custom_labels = [String(split(label, ": ", limit=2)[end]) for label in custom_labels]
    end

    hidden_matids = parameters.marker_plot_params.hidden_composition_matids
    hidden_colorbar_bins = nothing
    if !isnothing(hidden_matids)
        indices_to_remove = Set{Int}()
        for m in hidden_matids
            tick_pos = Float64(n_bin + 1 - m)
            idx = findfirst(t -> t == tick_pos, ticks)
            if !isnothing(idx)
                push!(indices_to_remove, idx)
            end
        end
        keep = [i for i in 1:length(ticks) if i ∉ indices_to_remove]
        custom_labels = custom_labels[keep]
        hidden_colorbar_bins = [n_bin + 1 - m for m in hidden_matids if 1 <= m <= n_bin]
        n_visible = length(keep)
        ticks = collect(1.0:Float64(n_visible))
    end

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
        hidden_colorbar_bins=hidden_colorbar_bins,
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
