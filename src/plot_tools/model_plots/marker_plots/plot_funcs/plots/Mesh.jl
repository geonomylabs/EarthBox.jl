module Mesh

import CairoMakie
import ....PlotParametersManager: PlotParameters
import ....GridPlotsManager.PlotScalarArraysManager: PlotScalarArrays
import ....PlotDtypes: AxesType

function plot_mesh(
    parameters::PlotParameters,
    scalar_arrays::PlotScalarArrays,
    axes::AxesType
)::Nothing
    _plot_mesh = parameters.marker_plot_params.plot_mesh
    mesh_line_width = parameters.marker_plot_params.mesh_line_width
    if _plot_mesh == 1
        println(">> Plotting mesh")
        for x in scalar_arrays.gridx
            CairoMakie.vlines!(
                axes, x;
                color=:black,
                linewidth=mesh_line_width
            )
        end
        for y in scalar_arrays.gridy
            CairoMakie.hlines!(
                axes, y;
                color=:black,
                linewidth=mesh_line_width
            )
        end
    end
    
    return nothing
end

end # module
