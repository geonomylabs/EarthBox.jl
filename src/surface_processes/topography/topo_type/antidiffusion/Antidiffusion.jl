module Antidiffusion

import EarthBox.ModelDataContainer: ModelData
import ....SedimentTransport: run_sediment_transport_model!
import ...AdvectionTools
import ...CarbonateAdvect

""" Advect topography vertically and horizontally.

Notes
-----
This module should be updated for the new downhill diffusion model and
cleaned up advection modules.

This should also be updated taking extrusion into account. Modify based on
the amount and location of extruded particles.
"""
function topography_advect!(
    model::ModelData,
    inside_flags::Vector{Int8}, # Not used
)::Nothing
    iuse_carb = model.carbonate.parameters.iuse_carb.value

    # Define advection velocity at topography grid points
    # (gridt[3,...], gridt[4,...])
    AdvectionTools.calc_advection_veloc_from_basic_grid!(
        model, model.topography.arrays.gridt.array)

    AdvectionTools.calc_advection_veloc_from_basic_grid!(
        model, model.carbonate.arrays.grid_carb.array)

    (
        iuse_downhill_diffusion
    ) = model.topography.parameters.model_options.iuse_downhill_diffusion.value

    if iuse_downhill_diffusion == 1
        AdvectionTools.advect_vertically!(
            model,
            use_stokes=true,
            use_uniform_erosion_and_deposition=false
        )
        AdvectionTools.advect_horizontally_with_antidiffusion!(model)
        run_sediment_transport_model!(model)
    else
        if iuse_carb == 0
            AdvectionTools.advect_vertically!(
                model,
                use_stokes=true,
                use_uniform_erosion_and_deposition=true
            )
        else
            CarbonateAdvect.advect_vertically_carbonate!(model)
        end
        AdvectionTools.advect_horizontally_with_antidiffusion!(model)
        if iuse_carb == 1
            # Carbonate array is advected horizontally since only the
            # topography grid was taken into account before. This is an odd
            # function since it does not use anti-diffusion whereas the
            # previous step for horizontal advection of topography did
            CarbonateAdvect.advect_carb!(model)
        end
    end

    return nothing
end

end # module 