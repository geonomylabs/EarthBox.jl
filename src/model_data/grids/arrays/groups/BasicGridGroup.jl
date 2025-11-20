module BasicGridGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.grids.arrays"
const GRP_NAME = "basic"

const ADATA = get_eb_arrays()

"""
    BasicGrid <: AbstractArrayGroup

Array group containing coordinates and spacing for basic grid.

# Fields
- `gridx_b::`[`GridArray1DState`](@ref): $(ADATA.gridx_b.description)
- `xstp_b::`[`GridArray1DState`](@ref): $(ADATA.xstp_b.description)
- `gridy_b::`[`GridArray1DState`](@ref): $(ADATA.gridy_b.description)
- `ystp_b::`[`GridArray1DState`](@ref): $(ADATA.ystp_b.description)

# Nested Dot Access
- `gridx_b = $(ROOT_NAME).$(GRP_NAME).gridx_b.array`
- `xstp_b = $(ROOT_NAME).$(GRP_NAME).xstp_b.array`
- `gridy_b = $(ROOT_NAME).$(GRP_NAME).gridy_b.array`
- `ystp_b = $(ROOT_NAME).$(GRP_NAME).ystp_b.array`

# Constructor
    BasicGrid(ynum::Int, xnum::Int)

Create a new BasicGrid array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `BasicGrid`: New BasicGrid array group with initialized arrays

"""
mutable struct BasicGrid <: AbstractArrayGroup
    gridx_b::GridArray1DState
    xstp_b::GridArray1DState
    gridy_b::GridArray1DState
    ystp_b::GridArray1DState
end

function BasicGrid(ynum::Int, xnum::Int)
    return BasicGrid(
        GridArray1DState(
            zeros(Float64, xnum),
            ADATA.gridx_b.name,
            ADATA.gridx_b.units,
            ADATA.gridx_b.description,
        ),
        GridArray1DState(
            zeros(Float64, xnum - 1),
            ADATA.xstp_b.name,
            ADATA.xstp_b.units,
            ADATA.xstp_b.description,
        ),
        GridArray1DState(
            zeros(Float64, ynum),
            ADATA.gridy_b.name,
            ADATA.gridy_b.units,
            ADATA.gridy_b.description,
        ),
        GridArray1DState(
            zeros(Float64, ynum - 1),
            ADATA.ystp_b.name,
            ADATA.ystp_b.units,
            ADATA.ystp_b.description,
        )
    )
end

end # module 