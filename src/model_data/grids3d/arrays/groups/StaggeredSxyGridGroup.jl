module StaggeredSxyGridGroup

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct StaggeredSxyGrid <: AbstractArrayGroup
    gridy_sxy::GridArray1DState
    ystp_sxy::GridArray1DState
    gridx_sxy::GridArray1DState
    xstp_sxy::GridArray1DState
    gridz_sxy::GridArray1DState
    zstp_sxy::GridArray1DState
end

function StaggeredSxyGrid(ynum::Int64, xnum::Int64, znum::Int64)
    return StaggeredSxyGrid(
        GridArray1DState(
            zeros(Float64, ynum), "gridy_sxy", "m",
            "y-locations of staggered deviatoric stress sxy grid nodes in x-y plane."
        ),
        GridArray1DState(
            zeros(Float64, ynum-1), "ystp_sxy", "m",
            "Width of staggered deviatoric stress sxy grid cells in y-direction."
        ),
        GridArray1DState(
            zeros(Float64, xnum), "gridx_sxy", "m",
            "x-locations of staggered deviatoric stress sxy grid nodes in x-y plane."
        ),
        GridArray1DState(
            zeros(Float64, xnum-1), "xstp_sxy", "m",
            "Width of staggered deviatoric stress sxy grid cells in x-direction."
        ),
        GridArray1DState(
            zeros(Float64, znum-1), "gridz_sxy", "m",
            "z-locations of staggered deviatoric stress sxy grid nodes in x-y plane."
        ),
        GridArray1DState(
            zeros(Float64, znum-2), "zstp_sxy", "m",
            "Width of staggered deviatoric stress sxy grid cells in z-direction."
        )
    )
end

end # module
