"""
    PlotSettingsManager

Module for managing general plot settings.

"""
module PlotSettingsManager

@Base.kwdef mutable struct PlotSettings
    plot_extension::String = ".png"
end

const PLOT_SETTINGS = PlotSettings()

end