module ImagePropsManager

import ...PlotDtypes: PlotDictType
import ...PlotDtypes: AbstractPlotParameterGroup

Base.@kwdef mutable struct ImageProps <: AbstractPlotParameterGroup
    aspect_ratio::Float64 = 1.0
    figure_dpi::Float64 = 150.0
    figsize::Tuple{Float64, Float64} = (5.0, 5.0)
    extension::String = ".png"
end

function ImageProps(plot_dict::PlotDictType)::ImageProps
    plot_params = plot_dict["general_parameters"]
    return ImageProps(
        aspect_ratio = get(plot_params, "aspect_ratio", 1.0),
        figure_dpi = get(plot_params, "figure_dpi", 150.0),
        figsize = get(plot_params, "figsize", (5.0, 5.0)),
        extension = get(plot_params, "extension", ".png")
    )
end

end # module