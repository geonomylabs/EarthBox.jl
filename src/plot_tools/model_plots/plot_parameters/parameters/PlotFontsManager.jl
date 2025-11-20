module PlotFontsManager

import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup

Base.@kwdef mutable struct PlotFonts <: AbstractPlotParameterGroup
    title_fontsize::Int = 12
    axis_title_fontsize::Int = 10
    axis_labels_fontsize::Int = 10
    axis_ticks_fontsize::Int = 10
    contour_label_fontsize::Int = 5
    legend_fontsize::Int = 6
    colorbar_ticks_fontsize::Int = 6
    colorbar_labels_fontsize::Int = 6
    text_box_font_size::Int = 6
    number_format::String = "%6.1f"
    rightside_up::Bool = true
end

function PlotFonts(plot_dict::PlotDictType)
    parameters = plot_dict["general_parameters"]
    return PlotFonts(
        title_fontsize=get(parameters, "title_fontsize", 12),
        axis_title_fontsize=get(parameters, "axis_title_fontsize", 10),
        axis_labels_fontsize=get(parameters, "axis_labels_fontsize", 10),
        axis_ticks_fontsize=get(parameters, "axis_ticks_fontsize", 10),
        contour_label_fontsize=get(parameters, "contour_label_fontsize", 5),
        legend_fontsize=get(parameters, "legend_fontsize", 6),
        colorbar_ticks_fontsize=get(parameters, "colorbar_ticks_fontsize", 6),
        colorbar_labels_fontsize=get(parameters, "colorbar_labels_fontsize", 6),
        text_box_font_size=get(parameters, "text_box_font_size", 6),
        number_format=get(parameters, "number_format", "%6.1f"),
        rightside_up=get(parameters, "rightside_up", true)
    )
end

function get_attribute_list(fonts::PlotFonts)::Vector{String}
    instance_attributes = fieldnames(typeof(fonts))
    return [string(attr) for attr in instance_attributes]
end

end # module
