module PlotPathsManager

import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup

Base.@kwdef mutable struct PlotPaths <: AbstractPlotParameterGroup
    mainpath::String = "None"
    outpath::String = "None"
end

function PlotPaths(plot_dict::PlotDictType)::PlotPaths
    plot_params = plot_dict["general_parameters"]
    return PlotPaths(
        mainpath = plot_params["mainpath"],
        outpath = plot_params["outpath"]
    )
end

end # module