module Initialize

import EarthBox.ModelDataContainer: ModelData

function initialize_stokes_global_loop_parameters!(model::ModelData)
    model.stokes_continuity.parameters.picard.iglobal.value = 0
    model.stokes_continuity.parameters.picard.iconverge.value = 0
    model.stokes_continuity.parameters.solution_norms.dvxy_rel_L2.value = 1e38
end

end # module 