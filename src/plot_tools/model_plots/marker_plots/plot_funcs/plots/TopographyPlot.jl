"""
    TopographyPlot

Module for plotting topography and base level in EarthBox visualization.
"""
module TopographyPlot

import ....PlotParametersManager: PlotParameters
import ...PlotTopoArraysManager: PlotTopoArrays
import ....PlotDtypes: AxesType
import CairoMakie

"""
    plot_topo(parameters::PlotParameters, topo_arrays::PlotTopoArrays, 
              axes::AxesType) -> Nothing

Plot topography on the given axes.

# Arguments
- `parameters::PlotParameters`: Plot parameters containing marker plot settings
- `topo_arrays::PlotTopoArrays`: Arrays containing topography data
- `axes::AxesType`: The axes object to plot on

# Returns
- `Nothing`: Function performs plotting operation
"""
function plot_topo(
    parameters::PlotParameters,
    topo_arrays::PlotTopoArrays,
    axes::AxesType
)::Nothing
    marker_plot_params = parameters.marker_plot_params
    if marker_plot_params.plot_topography == 1
        println(">> Plotting topography")
        topo_line_width = parameters.marker_plot_params.topo_line_width
        topo_line_color = parameters.marker_plot_params.topo_line_color
        CairoMakie.lines!(
            axes, topo_arrays.topox, topo_arrays.topoy;
            linewidth=topo_line_width, color=topo_line_color
        )
    end
    return nothing
end

"""
    plot_base_level(parameters::PlotParameters, topo_arrays::PlotTopoArrays, 
                   axes::AxesType) -> Nothing

Plot base level (sea level) line on the given axes.

# Arguments
- `parameters::PlotParameters`: Plot parameters containing marker plot settings
- `topo_arrays::PlotTopoArrays`: Arrays containing topography data
- `axes::AxesType`: The axes object to plot on

# Returns
- `Nothing`: Function performs plotting operation
"""
function plot_base_level(
    parameters::PlotParameters,
    topo_arrays::PlotTopoArrays,
    axes::AxesType
)::Nothing
    marker_plot_params = parameters.marker_plot_params
    y_sealevel = parameters.y_sealevel
    xmin = minimum(topo_arrays.topox)
    xmax = maximum(topo_arrays.topox)
    if marker_plot_params.plot_base_level == 1
        line_width = parameters.marker_plot_params.base_level_line_width
        line_color = parameters.marker_plot_params.base_level_line_color
        CairoMakie.lines!(
            axes, [xmin, xmax], [y_sealevel, y_sealevel];
            linewidth=line_width, color=line_color
        )
    end
    return nothing
end

end # module
