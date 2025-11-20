module MeanResidualsGroup

import EarthBox.Arrays.ArrayTypes.Array1DFloat: Array1DFloatState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct MeanResiduals <: AbstractArrayGroup
    resx_principle::Array1DFloatState
    resy_principle::Array1DFloatState
    resz_principle::Array1DFloatState
    resc_principle::Array1DFloatState
    resx_global::Array1DFloatState
    resy_global::Array1DFloatState
    resz_global::Array1DFloatState
    resc_global::Array1DFloatState
end

function MeanResiduals(nvcycles::Int)::MeanResiduals
    return MeanResiduals(
        Array1DFloatState(zeros(Float64, nvcycles), "resx_principle", "None", "Mean residual for x-Stokes equation for each multigrid iteration"),
        Array1DFloatState(zeros(Float64, nvcycles), "resy_principle", "None", "Mean residual for y-Stokes equation for each multigrid iteration"),
        Array1DFloatState(zeros(Float64, nvcycles), "resz_principle", "None", "Mean residual for z-Stokes equation for each multigrid iteration"),
        Array1DFloatState(zeros(Float64, nvcycles), "resc_principle", "None", "Mean residual for continuity equation for each multigrid iteration"),
        Array1DFloatState(zeros(Float64, nvcycles), "resx_global", "None", "Mean residual for x-Stokes equation for each multi-multigrid iteration"),
        Array1DFloatState(zeros(Float64, nvcycles), "resy_global", "None", "Mean residual for y-Stokes equation for each multi-multigrid iteration"),
        Array1DFloatState(zeros(Float64, nvcycles), "resz_global", "None", "Mean residual for z-Stokes equation for each multi-multigrid iteration"),
        Array1DFloatState(zeros(Float64, nvcycles), "resc_global", "None", "Mean residual for continuity equation for each multi-multigrid iteration")
    )
end

end