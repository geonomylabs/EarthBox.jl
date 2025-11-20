module StaggeredVxGridGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.grids.arrays"
const GRP_NAME = "staggered_vx"

const ADATA = get_eb_arrays()

"""
    StaggeredVxGrid <: AbstractArrayGroup

Array group containing coordinates and spacing for staggered Vx grid.

# Fields
- `gridy_vx::`[`GridArray1DState`](@ref): $(ADATA.gridy_vx.description)
- `ystp_vx::`[`GridArray1DState`](@ref): $(ADATA.ystp_vx.description)

# Nested Dot Access
- `gridy_vx = $(ROOT_NAME).$(GRP_NAME).gridy_vx.array`
- `ystp_vx = $(ROOT_NAME).$(GRP_NAME).ystp_vx.array`

# Constructor
    StaggeredVxGrid(ynum::Int)

Create a new StaggeredVxGrid array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction

# Returns
- `StaggeredVxGrid`: New StaggeredVxGrid array group with initialized arrays
"""
mutable struct StaggeredVxGrid <: AbstractArrayGroup
    gridy_vx::GridArray1DState
    ystp_vx::GridArray1DState
end

function StaggeredVxGrid(ynum::Int)
    return StaggeredVxGrid(
        GridArray1DState(
            zeros(Float64, ynum + 1),
            ADATA.gridy_vx.name,
            ADATA.gridy_vx.units,
            ADATA.gridy_vx.description,
        ),
        GridArray1DState(
            zeros(Float64, ynum),
            ADATA.ystp_vx.name,
            ADATA.ystp_vx.units,
            ADATA.ystp_vx.description,
        )
    )
end

end # module 