module PlasticStrainRateInterpolation

import EarthBox.ModelDataContainer: ModelData

"""
    interpolate_at_pressure_nodes!(model::ModelData)::Nothing

Interpolate plastic strain rate invariant at pressure nodes.

# Arguments
- `model::ModelData`: The model data structure containing all necessary arrays and parameters

# Updated Arrays
## Updated arrays from group `model.stokes_continuity.arrays.strain_rate_and_spin`:
- `eii_plastic_pressure.array::Matrix{Float64}` (ynum-1, xnum-1): 
    - Second invariant of plastic strain rate in 1/s.
"""
function interpolate_at_pressure_nodes!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    strain_rate_and_spin = model.stokes_continuity.arrays.strain_rate_and_spin
    eii_plastic_basic = strain_rate_and_spin.eii_plastic_basic.array
    eii_plastic_pressure = strain_rate_and_spin.eii_plastic_pressure.array

    Threads.@threads for j in 1:xnum-1
        for i in 1:ynum-1
            eii_upper_left = eii_plastic_basic[i,j]
            eii_upper_right = eii_plastic_basic[i,j+1]
            eii_lower_left = eii_plastic_basic[i+1,j]
            eii_lower_right = eii_plastic_basic[i+1,j+1]
            eii_plastic_pressure[i,j] = (
                eii_upper_left +
                eii_upper_right +
                eii_lower_left +
                eii_lower_right
            ) / 4.0
        end
    end
    return nothing
end

end # module 