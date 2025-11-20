module SetAxisTitles

import CairoMakie
import ....PlotParametersManager: PlotParameters

function set_xy_axis_titles!(
    axes_xy::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    set_font_size!(axes_xy, parameters)
    length_units = parameters.conversion.plot_units.length_units
    axes_xy.xlabel = "X ($length_units)"
    axes_xy.ylabel = "Y ($length_units)"
    return nothing
end

function set_heatflow_axis_titles!(
    axes_heatflow::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    set_font_size!(axes_heatflow, parameters)
    length_units = parameters.conversion.plot_units.length_units
    axes_heatflow.xlabel = "X ($length_units)"
    axes_heatflow.ylabel = "HeatFlow mW/m²"
    return nothing
end

function set_gravity_axis_titles!(
    axes_gravity::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    set_font_size!(axes_gravity, parameters)
    length_units = parameters.conversion.plot_units.length_units
    axes_gravity.xlabel = "X ($length_units)"
    axes_gravity.ylabel = "Gravity m/s²"
    return nothing
end

function set_font_size!(
    axes::CairoMakie.Axis, 
    parameters::PlotParameters
)::Nothing
    axes.titlesize = parameters.fonts.title_fontsize
    axes.xlabelsize = parameters.fonts.axis_title_fontsize
    axes.ylabelsize = parameters.fonts.axis_title_fontsize
    return nothing
end

end # module