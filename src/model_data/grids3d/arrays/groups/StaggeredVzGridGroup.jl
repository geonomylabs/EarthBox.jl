module StaggeredVzGridGroup

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct StaggeredVzGrid <: AbstractArrayGroup
    gridy_vz::GridArray1DState
    ystp_vz::GridArray1DState
    gridx_vz::GridArray1DState
    xstp_vz::GridArray1DState
    gridz_vz::GridArray1DState
    zstp_vz::GridArray1DState
end

function StaggeredVzGrid(ynum::Int64, xnum::Int64, znum::Int64)
    return StaggeredVzGrid(
        GridArray1DState(
            zeros(Float64, ynum+1), "gridy_vz", "m",
            "y-locations of staggered Vy grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, ynum), "ystp_vz", "m",
            "Width of staggered vy grid cells in y-direction."
        ),
        GridArray1DState(
            zeros(Float64, xnum + 1), "gridx_vz", "m",
            "x-locations of staggered Vy grid nodes including ghost nodes."
        ),
        GridArray1DState(
            zeros(Float64, xnum), "xstp_vz", "m",
            "Width of staggered vy grid cells in x-direction."
        ),
        GridArray1DState(
            zeros(Float64, znum), "gridz_vz", "m",
            "z-locations of staggered Vy grid nodes."
        ),
        GridArray1DState(
            zeros(Float64, znum-1), "zstp_vz", "m",
            "Width of staggered vy grid cells in z-direction."
        )
    )
end

end # module
