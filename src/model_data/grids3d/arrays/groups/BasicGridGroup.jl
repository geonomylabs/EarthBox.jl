module BasicGridGroup

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct BasicGrid <: AbstractArrayGroup
    gridy_b::GridArray1DState
    ystp_b::GridArray1DState
    gridx_b::GridArray1DState
    xstp_b::GridArray1DState
    gridz_b::GridArray1DState
    zstp_b::GridArray1DState
end

function BasicGrid(ynum::Int64, xnum::Int64, znum::Int64)
    return BasicGrid(
        GridArray1DState(
            zeros(Float64, ynum), "gridy_b", "m",
            "y-locations of basic grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, ynum - 1), "ystp_b", "m",
            "Width of cells in y-directions for basic grid."
        ),
        GridArray1DState(
            zeros(Float64, xnum), "gridx_b", "m",
            "x-locations of basic grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, xnum - 1), "xstp_b", "m",
            "Width of cells in x-directions for basic grid."
        ),
        GridArray1DState(
            zeros(Float64, znum), "gridz_b", "m",
            "z-locations of basic grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, znum - 1), "zstp_b", "m",
            "Width of cells in z-directions for basic grid."
        )
    )
end

end # module 