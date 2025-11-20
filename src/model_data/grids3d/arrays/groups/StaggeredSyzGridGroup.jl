module StaggeredSyzGridGroup

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct StaggeredSyzGrid <: AbstractArrayGroup
    gridy_syz::GridArray1DState
    ystp_syz::GridArray1DState
    gridx_syz::GridArray1DState
    xstp_syz::GridArray1DState
    gridz_syz::GridArray1DState
    zstp_syz::GridArray1DState
end

function StaggeredSyzGrid(ynum::Int64, xnum::Int64, znum::Int64)
    return StaggeredSyzGrid(
        GridArray1DState(
            zeros(Float64, ynum), "gridy_syz", "m",
            "y-locations of staggered deviatoric stress syz grid nodes in y-z plane."
        ),
        GridArray1DState(
            zeros(Float64, ynum-1), "ystp_syz", "m",
            "Width of staggered deviatoric stress syz grid cells in y-direction."
        ),
        GridArray1DState(
            zeros(Float64, xnum-1), "gridx_syz", "m",
            "x-locations of staggered deviatoric stress syz grid nodes in y-z plane."
        ),
        GridArray1DState(
            zeros(Float64, xnum-2), "xstp_syz", "m",
            "Width of staggered deviatoric stress syz grid cells in x-direction."
        ),
        GridArray1DState(
            zeros(Float64, znum), "gridz_syz", "m",
            "z-locations of staggered deviatoric stress syz grid nodes in y-z plane."
        ),
        GridArray1DState(
            zeros(Float64, znum-1), "zstp_syz", "m",
            "Width of staggered deviatoric stress syz grid cells in z-direction."
        )
    )
end

end # module
