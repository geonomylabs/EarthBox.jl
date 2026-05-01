""" Plot the lava thickness on the topography grid.
"""
module PlotLavaFlow

using CairoMakie
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
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "Y (m)", ylabel = "Thickness (m)")
    lines!(ax, topo_gridx, lava_thickness; color = :red)
    plot_name = string(
        "lava_thickness",
        "_drain_", idrainage_basin,
        "_time_", ntimestep,
        ".png"
    )
    save(plot_name, fig)
    return nothing
end

end # module 