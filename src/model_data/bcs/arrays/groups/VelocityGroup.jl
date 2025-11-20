"""
Module for combined velocity boundary condition arrays.

Provides data structures for specifying boundary conditions for all velocity
components together along each domain boundary.
"""
module VelocityGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.EarthBoxDtypes: AbstractArrayGroup
import EarthBox.Arrays.ArrayTypes.BcArrayFloat: BcArrayFloatState

const ROOT_NAME = "model.bcs.arrays"
const GRP_NAME = "velocity"

const ADATA = get_eb_arrays()

"""
    Velocity <: AbstractArrayGroup

Array group for combined velocity boundary conditions.

# Fields
- `btop::`[`BcArrayFloatState`](@ref): $(ADATA.btop.description)
- `bbottom::`[`BcArrayFloatState`](@ref): $(ADATA.bbottom.description)
- `bleft::`[`BcArrayFloatState`](@ref): $(ADATA.bleft.description)
- `bright::`[`BcArrayFloatState`](@ref): $(ADATA.bright.description)
- `bfront::`[`BcArrayFloatState`](@ref): $(ADATA.bfront.description)
- `bback::`[`BcArrayFloatState`](@ref): $(ADATA.bback.description)

# Nested Dot Access
- `btop = $(ROOT_NAME).$(GRP_NAME).btop.array`
- `bbottom = $(ROOT_NAME).$(GRP_NAME).bbottom.array`
- `bleft = $(ROOT_NAME).$(GRP_NAME).bleft.array`
- `bright = $(ROOT_NAME).$(GRP_NAME).bright.array`
- `bfront = $(ROOT_NAME).$(GRP_NAME).bfront.array`
- `bback = $(ROOT_NAME).$(GRP_NAME).bback.array`

# Constructor
    Velocity(ynum::Int64, xnum::Int64; znum::Int64=1)

Initializes combined velocity boundary condition arrays with zero values.

## Arguments
- `ynum::Int64`: Number of grid points in the vertical direction
- `xnum::Int64`: Number of grid points in the horizontal direction
- `znum::Int64=1`: Number of grid points in the z-direction (default: 1 for 2D)
"""
mutable struct Velocity <: AbstractArrayGroup
    btop::BcArrayFloatState
    bbottom::BcArrayFloatState
    bleft::BcArrayFloatState
    bright::BcArrayFloatState
    bfront::BcArrayFloatState
    bback::BcArrayFloatState
end

function Velocity(ynum::Int64, xnum::Int64; znum::Int64=1)::Velocity
    
    return Velocity(
        BcArrayFloatState(
            # Note that xnum +1 is used to account for the ghost nodes for the vy grid.
            # This allows a single multidimensional array to store the boundary conditions
            # for all three velocity components.
            zeros(Float64, xnum+1, 6), # array
            ADATA.btop.name, # name
            ADATA.btop.units, # units
            ADATA.btop.description
        ),
        BcArrayFloatState(
            zeros(Float64, xnum+1, 6), # array
            ADATA.bbottom.name, # name
            ADATA.bbottom.units, # units
            ADATA.bbottom.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 6), # array
            ADATA.bleft.name, # name
            ADATA.bleft.units, # units
            ADATA.bleft.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 6), # array
            ADATA.bright.name, # name
            ADATA.bright.units, # units
            ADATA.bright.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 6), # array
            ADATA.bfront.name, # name
            ADATA.bfront.units, # units
            ADATA.bfront.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 6), # array
            ADATA.bback.name, # name
            ADATA.bback.units, # units
            ADATA.bback.description
        )
    )
end

end # module 