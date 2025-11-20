module StaggeredVyGridGroup

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct StaggeredVyGrid <: AbstractArrayGroup
    gridy_vy::GridArray1DState
    ystp_vy::GridArray1DState
    gridx_vy::GridArray1DState
    xstp_vy::GridArray1DState
    gridz_vy::GridArray1DState
    zstp_vy::GridArray1DState
end

function StaggeredVyGrid(ynum::Int64, xnum::Int64, znum::Int64)
    return StaggeredVyGrid(
        GridArray1DState(
            zeros(Float64, ynum), "gridy_vy", "m",
            "y-locations of staggered Vy grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, ynum - 1), "ystp_vy", "m",
            "Width of staggered vy grid cells in y-direction."
        ),
        GridArray1DState(
            zeros(Float64, xnum + 1), "gridx_vy", "m",
            "x-locations of staggered Vy grid nodes including ghost nodes."
        ),
        GridArray1DState(
            zeros(Float64, xnum), "xstp_vy", "m",
            "Width of staggered vy grid cells in x-direction."
        ),
        GridArray1DState(
            zeros(Float64, znum + 1), "gridz_vy", "m",
            "z-locations of staggered Vy grid nodes including ghost nodes."
        ),
        GridArray1DState(
            zeros(Float64, znum), "zstp_vy", "m",
            "Width of staggered vy grid cells in z-direction."
        )
    )
end

end # module
