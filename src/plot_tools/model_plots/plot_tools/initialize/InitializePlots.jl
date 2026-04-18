module InitializePlots

import CairoMakie
import ..PlotUtils: get_figsize_pixels, normalize_height_ratios
import ..SetPlotAxes: set_xy_axes!, set_heatflow_axes!, set_gravity_axes!
import ...PlotParametersManager: PlotParameters

const margin_mm = 5

# Makie: (left, right, bottom, top). Extra top space for upper subplot y-tick labels.
const FIGURE_PADDING = (12, 12, 12, 28)

function initialize_xy_plot(
    parameters::PlotParameters
)::Tuple{CairoMakie.Figure, CairoMakie.Axis}
    fig, axes_xy = setup_xy_plot(parameters)
    set_xy_axes!(axes_xy, parameters)
    apply_layout_column_aspect_for_data_aspect!(fig, parameters)
    return fig, axes_xy
end

function apply_layout_column_aspect_for_data_aspect!(
    fig::CairoMakie.Figure,
    parameters::PlotParameters,
)::Nothing
    if !parameters.image.use_data_aspect
        return nothing
    end
    v = parameters.view
    dx = v.xmax_active - v.xmin_active
    dy = v.ymax_active - v.ymin_active
    if dx <= 0 || dy <= 0
        return nothing
    end
    ratio = dx / dy
    CairoMakie.colsize!(fig.layout, 1, CairoMakie.Aspect(1, Float32(ratio)))
    return nothing
end

function initialize_heatflow_composition_plot(
    parameters::PlotParameters
)::Tuple{CairoMakie.Figure, CairoMakie.Axis, CairoMakie.Axis}
    fig, axes_markers, axes_heatflow = setup_double_stacked_plot(parameters)
    set_xy_axes!(axes_markers, parameters)
    set_heatflow_axes!(axes_heatflow, parameters)
    return fig, axes_markers, axes_heatflow
end

function initialize_gravity_composition_plot(
    parameters::PlotParameters
)::Tuple{CairoMakie.Figure, CairoMakie.Axis, CairoMakie.Axis}
    fig, axes_markers, axes_gravity = setup_double_stacked_plot(parameters)
    set_xy_axes!(axes_markers, parameters)
    set_gravity_axes!(axes_gravity, parameters)
    return fig, axes_markers, axes_gravity
end

function initialize_heatflow_gravity_composition_plot(
    parameters::PlotParameters
)::Tuple{CairoMakie.Figure, CairoMakie.Axis, CairoMakie.Axis, CairoMakie.Axis}
    fig, axes_markers, axes_heatflow, axes_gravity = setup_triple_stacked_plot(parameters)
    set_xy_axes!(axes_markers, parameters)
    set_heatflow_axes!(axes_heatflow, parameters)
    set_gravity_axes!(axes_gravity, parameters)
    return fig, axes_markers, axes_heatflow, axes_gravity
end

function setup_xy_plot(parameters::PlotParameters)::Tuple{CairoMakie.Figure, CairoMakie.Axis}
    fig = CairoMakie.Figure(
        size=get_figsize_pixels(parameters),
        figure_padding=FIGURE_PADDING,
        backgroundcolor=:white
        )
    ax_aspect = parameters.image.use_data_aspect ? CairoMakie.DataAspect() : nothing
    axes_xy = CairoMakie.Axis(fig[1, 1], yreversed=true, aspect=ax_aspect)
    return fig, axes_xy
end

""" Setup a plot with two stacked subplots.

The plot on the top is for a curve plot (e.g. heat flow or gravity) and the plot
on the bottom is for a scatter plot of markers.
"""
function setup_double_stacked_plot(
    parameters::PlotParameters
)::Tuple{CairoMakie.Figure, CairoMakie.Axis, CairoMakie.Axis}
    height_ratios = parameters.marker_plot_params.height_ratios
    if isnothing(height_ratios)
        height_ratios = [0.25, 0.75]
    elseif length(height_ratios) < 2
        height_ratios = [0.25, 0.75]
    elseif length(height_ratios) > 2
        height_ratios = [0.25, 0.75]
    end
    height_ratios = normalize_height_ratios(height_ratios)
    fig = CairoMakie.Figure(size=get_figsize_pixels(parameters), figure_padding=FIGURE_PADDING)
    axes_scatter = CairoMakie.Axis(fig[2, 1], yreversed=true)
    axes_curve = CairoMakie.Axis(fig[1, 1], yreversed=false)
    fig.layout.rowsizes[1] = CairoMakie.Relative(height_ratios[1])
    fig.layout.rowsizes[2] = CairoMakie.Relative(height_ratios[2])
    return fig, axes_scatter, axes_curve
end

""" Setup a plot with three stacked subplots.

The plot on the top and middle are for a curve plots (e.g. heat flow or gravity)
and the plot on the bottom is for a scatter plot of markers.
"""
function setup_triple_stacked_plot(
    parameters::PlotParameters
)::Tuple{CairoMakie.Figure, CairoMakie.Axis, CairoMakie.Axis, CairoMakie.Axis}
    height_ratios = parameters.marker_plot_params.height_ratios
    if isnothing(height_ratios)
        height_ratios = [0.25, 0.25, 0.5]
    elseif length(height_ratios) < 3
        height_ratios = [0.25, 0.25, 0.5]
    elseif length(height_ratios) > 3
        height_ratios = [0.25, 0.25, 0.5]
    end
    height_ratios = normalize_height_ratios(height_ratios)
    fig = CairoMakie.Figure(size=get_figsize_pixels(parameters), figure_padding=FIGURE_PADDING)
    axes_scatter = CairoMakie.Axis(fig[3, 1], yreversed=true)
    axes_curve1 = CairoMakie.Axis(fig[2, 1], yreversed=false)
    axes_curve2 = CairoMakie.Axis(fig[1, 1], yreversed=false)
    fig.layout.rowsizes[1] = CairoMakie.Relative(height_ratios[1])
    fig.layout.rowsizes[2] = CairoMakie.Relative(height_ratios[2])
    fig.layout.rowsizes[3] = CairoMakie.Relative(height_ratios[3])
    return fig, axes_scatter, axes_curve1, axes_curve2
end

end