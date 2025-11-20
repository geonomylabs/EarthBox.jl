module HarmonicAverage

import EarthBox.ModelDataContainer: ModelData
import EarthBox: MathTools
import ..CheckPlasticYield: plastic_yielding

"""
    normal_viscosity_from_harmonic_avg_shear_viscosity!(model::ModelData)::Nothing

Calculate normal viscosity with a harmonic average of shear viscosity.

# Arguments
- `model::ModelData`: The model data structure containing all necessary arrays and parameters

# Updated Arrays
## Updated arrays from group `model.stokes_continuity.arrays.viscosity`:
- `etan1.array::Matrix{Float64}` (ynum-1, xnum-1): 
    - Normal viscosity (Pa.s) interpolated from markers to pressure grid.
"""
function normal_viscosity_from_harmonic_avg_shear_viscosity!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    etan1 = model.stokes_continuity.arrays.viscosity.etan1.array
    etas1 = model.stokes_continuity.arrays.viscosity.etas1.array
    plastics = model.stokes_continuity.arrays.plastic_def.plastics.array
    plasticn = model.stokes_continuity.arrays.plastic_def.plasticn.array
    plastic_yield = model.stokes_continuity.arrays.plastic_def.plastic_yield.array
    icheck = 0
    plas_tol = 0.0

    Threads.@threads for i in 1:ynum-1
        for j in 1:xnum-1
            yielding = plastic_yielding(
                i, j, plasticn, plastics, plastic_yield, plas_tol
            )
            if yielding
                eta_update, icheck = MathTools.harmonic4(
                    etas1[i,j],
                    etas1[i,j+1],
                    etas1[i+1,j],
                    etas1[i+1,j+1],
                    0.0, 0.0, 0.0, 0.0
                )
                if icheck != 1
                    etan1[i,j] = eta_update
                end
            end
        end
    end
    return nothing
end

end # module 