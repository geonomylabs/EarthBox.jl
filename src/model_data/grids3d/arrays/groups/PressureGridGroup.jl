module PressureGridGroup

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct PressureGrid <: AbstractArrayGroup
    gridy_pr::GridArray1DState
    ystp_pr::GridArray1DState
    gridx_pr::GridArray1DState
    xstp_pr::GridArray1DState
    gridz_pr::GridArray1DState
    zstp_pr::GridArray1DState
end

function PressureGrid(ynum::Int64, xnum::Int64, znum::Int64)
    return PressureGrid(
        GridArray1DState(
            zeros(Float64, ynum - 1), "gridy_pr", "m", 
            "y-locations of pressure grid nodes."
            ),
        GridArray1DState(
            zeros(Float64, ynum - 2), "ystp_pr", "m",
            "Width of pressure cells in y-direction."
        ),
        GridArray1DState(
            zeros(Float64, xnum - 1), "gridx_pr", "m",
            "x-locations of pressure grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, xnum - 2), "xstp_pr", "m",
            "Width of pressure cells in x-direction."
        ),
        GridArray1DState(
            zeros(Float64, znum - 1), "gridz_pr", "m",
            "z-locations of pressure grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, znum - 2), "zstp_pr", "m",
            "Width of pressure cells in z-direction."
        )
    )
end

end # module 