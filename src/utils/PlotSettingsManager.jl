"""
    PlotSettingsManager

Module for managing general plot settings.

"""
module PlotSettingsManager

function __init__()
    get!(ENV, "GKSwstype", "100")
end

@Base.kwdef mutable struct PlotSettings
    plot_extension::String = ".png"
end

const PLOT_SETTINGS = PlotSettings()

end