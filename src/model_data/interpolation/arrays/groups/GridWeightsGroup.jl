module GridWeightsGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.interpolation.arrays"
const GRP_NAME = "grid_weights"

const ADATA = get_eb_arrays()

"""
    GridWeights <: AbstractArrayGroup

Array group containing summed interpolation weight arrays on grid nodes.

# Fields
- `wtnodes::`[`ScalarArray2DState`](@ref): $(ADATA.wtnodes.description)
- `wtetas::`[`ScalarArray2DState`](@ref): $(ADATA.wtetas.description)
- `wtetan::`[`ScalarArray2DState`](@ref): $(ADATA.wtetan.description)
- `wtnodes_vy::`[`ScalarArray2DState`](@ref): $(ADATA.wtnodes_vy.description)

# Nested Dot Access
- `wtnodes = $(ROOT_NAME).$(GRP_NAME).wtnodes.array`
- `wtetas = $(ROOT_NAME).$(GRP_NAME).wtetas.array`
- `wtetan = $(ROOT_NAME).$(GRP_NAME).wtetan.array`
- `wtnodes_vy = $(ROOT_NAME).$(GRP_NAME).wtnodes_vy.array`

# Constructor
    GridWeights(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `GridWeights`: New GridWeights array group with initialized arrays

"""
mutable struct GridWeights <: AbstractArrayGroup
    wtnodes::ScalarArray2DState
    wtetas::ScalarArray2DState
    wtetan::ScalarArray2DState
    wtnodes_vy::ScalarArray2DState
end

function GridWeights(ynum::Int, xnum::Int)::GridWeights
    return GridWeights(
        ScalarArray2DState(
            ynum,                                 # ynum
            xnum,                                 # xnum
            ADATA.wtnodes.name,                  # name
            ADATA.wtnodes.units,                 # units
            ADATA.wtnodes.grid_type,             # grid_type
            ADATA.wtnodes.description            # description
        ),
        ScalarArray2DState(
            ynum,                                 # ynum
            xnum,                                 # xnum
            ADATA.wtetas.name,                   # name
            ADATA.wtetas.units,                  # units
            ADATA.wtetas.grid_type,              # grid_type
            ADATA.wtetas.description             # description
        ),
        ScalarArray2DState(
            ynum,                                 # ynum
            xnum,                                 # xnum
            ADATA.wtetan.name,                   # name
            ADATA.wtetan.units,                  # units
            ADATA.wtetan.grid_type,              # grid_type
            ADATA.wtetan.description             # description
        ),
        ScalarArray2DState(
            ynum,                                 # ynum
            xnum,                                 # xnum
            ADATA.wtnodes_vy.name,               # name
            ADATA.wtnodes_vy.units,              # units
            ADATA.wtnodes_vy.grid_type,          # grid_type
            ADATA.wtnodes_vy.description         # description
        )
    )
end

end # module 