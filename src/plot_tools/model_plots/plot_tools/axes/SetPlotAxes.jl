module SetPlotAxes

include("core/SetAxisTitles.jl")
include("core/SetPlotView.jl")
include("core/SetPlotTicks.jl")

import CairoMakie
import .SetAxisTitles: set_xy_axis_titles!, set_heatflow_axis_titles!, set_gravity_axis_titles!
import .SetPlotView: set_xy_view!, set_heatflow_view!, set_gravity_view!
import .SetPlotTicks: set_xy_ticks!, set_heatflow_ticks!, set_gravity_ticks!
import ...PlotParametersManager: PlotParameters

function set_xy_axes!(
    axes_xy::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    set_xy_view!(axes_xy, parameters)
    set_xy_ticks!(axes_xy, parameters)
    set_xy_axis_titles!(axes_xy, parameters)
    return nothing
end

function set_heatflow_axes!(
    axes_heatflow::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    set_heatflow_view!(axes_heatflow, parameters)
    set_heatflow_ticks!(axes_heatflow, parameters)
    set_heatflow_axis_titles!(axes_heatflow, parameters)
    return nothing
end

function set_gravity_axes!(
    axes_gravity::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    set_gravity_view!(axes_gravity, parameters)
    set_gravity_ticks!(axes_gravity, parameters)
    set_gravity_axis_titles!(axes_gravity, parameters)
    return nothing
end

end # module