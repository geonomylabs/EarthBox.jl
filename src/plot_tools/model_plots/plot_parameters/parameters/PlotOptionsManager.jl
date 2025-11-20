module PlotOptionsManager

import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup

Base.@kwdef mutable struct PlotOptions <: AbstractPlotParameterGroup
    iplot::Int = 0
    grid_plot_type::String = "None"
    linewidth::Float64 = 0.0
    edgecolor::String = "None"
    show_nodes::Int = 0
end

function PlotOptions(plot_dict::PlotDictType)::PlotOptions
    plot_params = plot_dict["general_parameters"]
    return PlotOptions(
        iplot = get(plot_params, "iplot", 0),
        grid_plot_type = get(plot_params, "grid_plot_type", "None"),
        linewidth = get(plot_params, "linewidth", 0.0),
        edgecolor = get(plot_params, "edgecolor", "None"),
        show_nodes = get(plot_params, "show_nodes", 0)
    )
end

end # module