module PlotTicksManager

import ..PlotViewManager: PlotView
import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup

Base.@kwdef mutable struct PlotTicks <: AbstractPlotParameterGroup
    xspacing::Float64 = 0.0
    yspacing::Float64 = 0.0
    xticks::Vector{Float64} = [0.0]
    yticks::Vector{Float64} = [0.0]
end

function PlotTicks(plot_dict::PlotDictType)
    params = plot_dict["general_parameters"]
    default_tuple = (0.0, 0.0)
    xspacing = get(params, "xyspacing", default_tuple)[1]
    yspacing = get(params, "xyspacing", default_tuple)[2]
    return PlotTicks(
        xspacing=xspacing,
        yspacing=yspacing,
        xticks=[0.0],
        yticks=[0.0],
    )
end

function update_ticks!(ticks::PlotTicks, view::PlotView)::Nothing
    ticks.xticks = calculate_axis_ticks(
        view.xmin_active, view.xmax_active, ticks.xspacing)
    ticks.yticks = calculate_axis_ticks(
        view.ymin_active, view.ymax_active, ticks.yspacing)
    return nothing
end

function calculate_axis_ticks(
    axis_value_min::Float64,
    axis_value_max::Float64,
    tick_width::Float64
)::Vector{Float64}
    axis_size = axis_value_max - axis_value_min
    nticks = floor(Int, axis_size/tick_width) + 1
    ticks = Float64[]
    for i in 1:nticks
        tick_value = axis_value_min + tick_width * (i - 1)
        push!(ticks, tick_value)
    end
    return ticks
end

end # module