module StaggeredVyGridGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.grids.arrays"
const GRP_NAME = "staggered_vy"

const ADATA = get_eb_arrays()

"""
    StaggeredVyGrid <: AbstractArrayGroup

Array group containing coordinates and spacing for staggered Vy grid.

# Fields
- `gridx_vy::`[`GridArray1DState`](@ref): $(ADATA.gridx_vy.description)
- `xstp_vy::`[`GridArray1DState`](@ref): $(ADATA.xstp_vy.description)

# Nested Dot Access
- `gridx_vy = $(ROOT_NAME).$(GRP_NAME).gridx_vy.array`
- `xstp_vy = $(ROOT_NAME).$(GRP_NAME).xstp_vy.array`

# Constructor
    StaggeredVyGrid(xnum::Int)

Create a new StaggeredVyGrid array group with the given grid dimensions.

# Arguments
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `StaggeredVyGrid`: New StaggeredVyGrid array group with initialized arrays
"""
mutable struct StaggeredVyGrid <: AbstractArrayGroup
    gridx_vy::GridArray1DState
    xstp_vy::GridArray1DState
end

function StaggeredVyGrid(xnum::Int)
    return StaggeredVyGrid(
        GridArray1DState(
            zeros(Float64, xnum + 1),
            ADATA.gridx_vy.name,
            ADATA.gridx_vy.units,
            ADATA.gridx_vy.description,
        ),
        GridArray1DState(
            zeros(Float64, xnum),
            ADATA.xstp_vy.name,
            ADATA.xstp_vy.units,
            ADATA.xstp_vy.description,
        )
    )
end

end # module
