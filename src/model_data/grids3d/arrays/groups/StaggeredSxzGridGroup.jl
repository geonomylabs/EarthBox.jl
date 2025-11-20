module StaggeredSxzGridGroup

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct StaggeredSxzGrid <: AbstractArrayGroup
    gridy_sxz::GridArray1DState
    ystp_sxz::GridArray1DState
    gridx_sxz::GridArray1DState
    xstp_sxz::GridArray1DState
    gridz_sxz::GridArray1DState
    zstp_sxz::GridArray1DState
end

function StaggeredSxzGrid(ynum::Int64, xnum::Int64, znum::Int64)
    return StaggeredSxzGrid(
        GridArray1DState(
            zeros(Float64, ynum-1), "gridy_sxz", "m",
            "y-locations of staggered deviatoric stress sxz grid nodes in x-z plane."
        ),
        GridArray1DState(
            zeros(Float64, ynum-2), "ystp_sxz", "m",
            "Width of staggered deviatoric stress sxz grid cells in y-direction."
        ),
        GridArray1DState(
            zeros(Float64, xnum), "gridx_sxz", "m",
            "x-locations of staggered deviatoric stress sxz grid nodes in x-z plane."
        ),
        GridArray1DState(
            zeros(Float64, xnum-1), "xstp_sxz", "m",
            "Width of staggered deviatoric stress sxz grid cells in x-direction."
        ),
        GridArray1DState(
            zeros(Float64, znum), "gridz_sxz", "m",
            "z-locations of staggered deviatoric stress sxz grid nodes in x-z plane."
        ),
        GridArray1DState(
            zeros(Float64, znum-1), "zstp_sxz", "m",
            "Width of staggered deviatoric stress sxz grid cells in z-direction."
        )
    )
end

end # module
