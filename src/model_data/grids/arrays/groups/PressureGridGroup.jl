module PressureGridGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.grids.arrays"
const GRP_NAME = "pressure"

const ADATA = get_eb_arrays()

"""
    PressureGrid <: AbstractArrayGroup

Array group containing coordinates and spacing for pressure grid.

# Fields
- `gridy_pr::`[`GridArray1DState`](@ref): $(ADATA.gridy_pr.description)
- `ystp_pr::`[`GridArray1DState`](@ref): $(ADATA.ystp_pr.description)
- `gridx_pr::`[`GridArray1DState`](@ref): $(ADATA.gridx_pr.description)
- `xstp_pr::`[`GridArray1DState`](@ref): $(ADATA.xstp_pr.description)

# Nested Dot Access
- `gridy_pr = $(ROOT_NAME).$(GRP_NAME).gridy_pr.array`
- `ystp_pr = $(ROOT_NAME).$(GRP_NAME).ystp_pr.array`
- `gridx_pr = $(ROOT_NAME).$(GRP_NAME).gridx_pr.array`
- `xstp_pr = $(ROOT_NAME).$(GRP_NAME).xstp_pr.array`

# Constructor
    PressureGrid(ynum::Int, xnum::Int)

Create a new PressureGrid array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `PressureGrid`: New PressureGrid array group with initialized arrays
"""
mutable struct PressureGrid <: AbstractArrayGroup
    gridy_pr::GridArray1DState
    ystp_pr::GridArray1DState
    gridx_pr::GridArray1DState
    xstp_pr::GridArray1DState
end

function PressureGrid(ynum::Int, xnum::Int)
    return PressureGrid(
        GridArray1DState(
            zeros(Float64, ynum - 1),
            ADATA.gridy_pr.name,
            ADATA.gridy_pr.units,
            ADATA.gridy_pr.description,
        ),
        GridArray1DState(
            zeros(Float64, ynum - 2),
            ADATA.ystp_pr.name,
            ADATA.ystp_pr.units,
            ADATA.ystp_pr.description,
        ),
        GridArray1DState(
            zeros(Float64, xnum - 1),
            ADATA.gridx_pr.name,
            ADATA.gridx_pr.units,
            ADATA.gridx_pr.description,
        ),
        GridArray1DState(
            zeros(Float64, xnum - 2),
            ADATA.xstp_pr.name,
            ADATA.xstp_pr.units,
            ADATA.xstp_pr.description,
        )
    )
end

end # module 