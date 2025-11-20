module SetPlotView

import CairoMakie
import ....PlotParametersManager: PlotParameters

function set_xy_view!(
    axes_xy::CairoMakie.Axis,
    parameters::PlotParameters;
)::Nothing
    view = parameters.view
    axes_xy.limits = (view.xmin_active, view.xmax_active, view.ymin_active, view.ymax_active)
    return nothing
end

function set_heatflow_view!(
    axes_heatflow::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    view = parameters.view
    min_hf = parameters.marker_plot_params.heatflow_min
    max_hf = parameters.marker_plot_params.heatflow_max
    axes_heatflow.limits = (view.xmin_active, view.xmax_active, min_hf, max_hf)
    return nothing
end

function set_gravity_view!(
    axes_gravity::CairoMakie.Axis,
    parameters::PlotParameters
)::Nothing
    view = parameters.view
    min_grav = parameters.marker_plot_params.gravity_min
    max_grav = parameters.marker_plot_params.gravity_max
    axes_gravity.limits = (view.xmin_active, view.xmax_active, min_grav, max_grav)
    return nothing
end

end # module