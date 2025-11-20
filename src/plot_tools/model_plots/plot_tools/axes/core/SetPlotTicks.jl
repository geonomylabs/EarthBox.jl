module SetPlotTicks

import CairoMakie
import ....PlotParametersManager: PlotParameters

function set_xy_ticks!(
    axes_xy::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    axes_xy.xticks = parameters.ticks.xticks
    axes_xy.yticks = parameters.ticks.yticks
    set_tick_fontsize!(axes_xy, parameters)
    return nothing
end

function set_heatflow_ticks!(
    axes_heatflow::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    axes_heatflow.xticks = parameters.ticks.xticks

    min_hf = parameters.marker_plot_params.heatflow_min
    max_hf = parameters.marker_plot_params.heatflow_max
    hf_spacing = parameters.marker_plot_params.heatflow_spacing
    hf_ticks = calculate_ticks(min_hf, max_hf, hf_spacing)

    axes_heatflow.yticks = hf_ticks
    set_tick_fontsize!(axes_heatflow, parameters)
    return nothing
end

function set_gravity_ticks!(
    axes_gravity::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    axes_gravity.xticks = parameters.ticks.xticks

    min_grav = parameters.marker_plot_params.gravity_min
    max_grav = parameters.marker_plot_params.gravity_max
    grav_spacing = parameters.marker_plot_params.gravity_spacing
    grav_ticks = calculate_ticks(min_grav, max_grav, grav_spacing)

    axes_gravity.yticks = grav_ticks
    set_tick_fontsize!(axes_gravity, parameters)
    return nothing
end

function set_tick_fontsize!(
    axes::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    tick_fontsize = parameters.fonts.axis_ticks_fontsize
    axes.xticklabelsize = tick_fontsize
    axes.yticklabelsize = tick_fontsize
    return nothing
end

function calculate_ticks(
    min_val::Float64,
    max_val::Float64,
    spacing::Float64
)::Vector{Float64}
    return [i for i in floor(Int, min_val):floor(Int, spacing):floor(Int, max_val)]
end

end # module