module RungeKutta

include("core/TopoFromRungeKutta.jl")

import EarthBox.ModelDataContainer: ModelData
import ...AdvectionTools: calc_advection_veloc_from_basic_grid!
import ...AdvectionTools: advect_vertically!
import .TopoFromRungeKutta: calc_topo_from_runge_kutta!

""" Advect topography vertically and horizontally.

The topography marker-chain is updated and old sediment markers are advected due to 
compaction of newly deposited sediment. Markers are not transformed in this function to 
account for changes in topography. These transformations are managed by the 
surface_processes module after the topography is updated.
"""
function topography_advect!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    calc_advection_veloc_from_basic_grid!(model, model.topography.arrays.gridt.array)
    runge_kutta_order_max = 4
    calc_topo_from_runge_kutta!(model, inside_flags, runge_kutta_order_max)
    (
        iuse_downhill_diffusion
    ) = model.topography.parameters.model_options.iuse_downhill_diffusion.value

    if iuse_downhill_diffusion == 0
        advect_vertically!(
            model,
            use_stokes=false,
            use_uniform_erosion_and_deposition=true
        )
    end
    
    return nothing
end

end # module 