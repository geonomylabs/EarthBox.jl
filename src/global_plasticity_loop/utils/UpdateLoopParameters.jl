module UpdateLoopParameters

import EarthBox.ModelDataContainer: ModelData

""" Update global loop parameters in main model interface.

# Updated parameters from group `model.stokes_continuity.parameters.picard`
- `iconverge.value::Int64`: Integer flag indicating convergence.
- `iglobal.value::Int64`: Global Picard loop iteration counter.
"""
function update!(model::ModelData)
    stokes_solution_norms = model.stokes_continuity.parameters.solution_norms
    tol_global = model.stokes_continuity.parameters.picard.tolerance_picard.value
    
    if stokes_solution_norms.dvxy_rel_L2.value <= tol_global
        model.stokes_continuity.parameters.picard.iconverge.value = 1
    end
    model.stokes_continuity.parameters.picard.iglobal.value += 1
end

end # module 