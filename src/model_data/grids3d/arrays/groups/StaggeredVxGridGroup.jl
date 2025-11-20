module StaggeredVxGridGroup

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct StaggeredVxGrid <: AbstractArrayGroup
    gridy_vx::GridArray1DState
    ystp_vx::GridArray1DState
    gridx_vx::GridArray1DState
    xstp_vx::GridArray1DState
    gridz_vx::GridArray1DState
    zstp_vx::GridArray1DState
end

function StaggeredVxGrid(ynum::Int64, xnum::Int64, znum::Int64)
    return StaggeredVxGrid(
        GridArray1DState(
            zeros(Float64, ynum + 1), "gridy_vx", "m",
            "y-locations of staggered Vx grid nodes including ghost nodes."
        ),
        GridArray1DState(
            zeros(Float64, ynum), "ystp_vx", "m",
            "Width of staggered vx grid cells in y-direction."
        ),
        GridArray1DState(
            zeros(Float64, xnum), "gridx_vx", "m",
            "x-locations of staggered Vx grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, xnum - 1), "xstp_vx", "m",
            "Width of staggered vx grid cells in x-direction."
        ),
        GridArray1DState(
            zeros(Float64, znum + 1), "gridz_vx", "m",
            "z-locations of staggered Vx grid nodes including ghost nodes."
        ),
        GridArray1DState(
            zeros(Float64, znum), "zstp_vx", "m",
            "Width of staggered vx grid cells in z-direction."
        )
    )
end

end # module 