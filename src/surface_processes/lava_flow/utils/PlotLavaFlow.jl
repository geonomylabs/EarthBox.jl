""" Plot the lava thickness on the topography grid.
"""
module PlotLavaFlow

import Plots
import EarthBox.ModelDataContainer: ModelData

""" Plot the lava thickness on the topography grid.
"""
function plot_lava_thickness(
    model::ModelData,
    topo_gridx::Vector{Float64},
    lava_thickness::Vector{Float64},
    idrainage_basin::Int
)::Nothing
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    p = Plots.plot(
        topo_gridx,
        lava_thickness,
        label = "",
        linecolor = :red,
        xlabel = "Y (m)",
        ylabel = "Thickness (m)"
    )
    plot_name = string(
        "lava_thickness",
        "_drain_", idrainage_basin,
        "_time_", ntimestep,
        ".png"
    )
    Plots.savefig(p, plot_name)
    return nothing
end

end # module 