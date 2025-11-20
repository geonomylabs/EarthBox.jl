module RemainingStress

import EarthBox.ModelDataContainer: ModelData

"""
    calculate_remaining_stress_change!(model::ModelData)

Calculate remaining stress change by subtracting relaxed nodal-marker stress 
difference from viscoelastic stress change:
    dsxy = dsxy - dsxyn
    dsxx = dsxx - dsxxn

# Updated arrays from group model.stokes_continuity.arrays.stress_change

- dsxy: Final (remaining) shear stress change (Pa) on basic grid.
- dsxx: Final (remaining) normal stress change (Pa) on pressure grid.
"""
function calculate_remaining_stress_change!(model::ModelData)::Nothing
    model.stokes_continuity.arrays.stress_change.dsxy.array .-= 
        model.stokes_continuity.arrays.stress_change.dsxyn.array

    model.stokes_continuity.arrays.stress_change.dsxx.array .-= 
        model.stokes_continuity.arrays.stress_change.dsxxn.array

    return nothing
end

end # module 