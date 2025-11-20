"""
Module for component-wise velocity boundary condition arrays.

Provides data structures for specifying boundary conditions for each velocity
component (x, y, z) separately along all domain boundaries.
"""
module VelCompGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.EarthBoxDtypes: AbstractArrayGroup
import EarthBox.Arrays.ArrayTypes.BcArrayFloat: BcArrayFloatState

const ROOT_NAME = "model.bcs.arrays"
const GRP_NAME = "vel_comp"

const ADATA = get_eb_arrays()

"""
    VelComp <: AbstractArrayGroup

Array group for component-wise velocity boundary conditions.

# Fields
- `btopx::`[`BcArrayFloatState`](@ref): $(ADATA.btopx.description)
- `btopy::`[`BcArrayFloatState`](@ref): $(ADATA.btopy.description)
- `btopz::`[`BcArrayFloatState`](@ref): $(ADATA.btopz.description)
- `bbottomx::`[`BcArrayFloatState`](@ref): $(ADATA.bbottomx.description)
- `bbottomy::`[`BcArrayFloatState`](@ref): $(ADATA.bbottomy.description)
- `bbottomz::`[`BcArrayFloatState`](@ref): $(ADATA.bbottomz.description)
- `bleftx::`[`BcArrayFloatState`](@ref): $(ADATA.bleftx.description)
- `blefty::`[`BcArrayFloatState`](@ref): $(ADATA.blefty.description)
- `bleftz::`[`BcArrayFloatState`](@ref): $(ADATA.bleftz.description)
- `brightx::`[`BcArrayFloatState`](@ref): $(ADATA.brightx.description)
- `brighty::`[`BcArrayFloatState`](@ref): $(ADATA.brighty.description)
- `brightz::`[`BcArrayFloatState`](@ref): $(ADATA.brightz.description)
- `bfrontx::`[`BcArrayFloatState`](@ref): $(ADATA.bfrontx.description)
- `bfronty::`[`BcArrayFloatState`](@ref): $(ADATA.bfronty.description)
- `bfrontz::`[`BcArrayFloatState`](@ref): $(ADATA.bfrontz.description)
- `bbackx::`[`BcArrayFloatState`](@ref): $(ADATA.bbackx.description)
- `bbacky::`[`BcArrayFloatState`](@ref): $(ADATA.bbacky.description)
- `bbackz::`[`BcArrayFloatState`](@ref): $(ADATA.bbackz.description)

# Nested Dot Access
- `btopx = $(ROOT_NAME).$(GRP_NAME).btopx.array`
- `btopy = $(ROOT_NAME).$(GRP_NAME).btopy.array`
- `btopz = $(ROOT_NAME).$(GRP_NAME).btopz.array`
- `bbottomx = $(ROOT_NAME).$(GRP_NAME).bbottomx.array`
- `bbottomy = $(ROOT_NAME).$(GRP_NAME).bbottomy.array`
- `bbottomz = $(ROOT_NAME).$(GRP_NAME).bbottomz.array`
- `bleftx = $(ROOT_NAME).$(GRP_NAME).bleftx.array`
- `blefty = $(ROOT_NAME).$(GRP_NAME).blefty.array`
- `bleftz = $(ROOT_NAME).$(GRP_NAME).bleftz.array`
- `brightx = $(ROOT_NAME).$(GRP_NAME).brightx.array`
- `brighty = $(ROOT_NAME).$(GRP_NAME).brighty.array`
- `brightz = $(ROOT_NAME).$(GRP_NAME).brightz.array`
- `bfrontx = $(ROOT_NAME).$(GRP_NAME).bfrontx.array`
- `bfronty = $(ROOT_NAME).$(GRP_NAME).bfronty.array`
- `bfrontz = $(ROOT_NAME).$(GRP_NAME).bfrontz.array`
- `bbackx = $(ROOT_NAME).$(GRP_NAME).bbackx.array`
- `bbacky = $(ROOT_NAME).$(GRP_NAME).bbacky.array`
- `bbackz = $(ROOT_NAME).$(GRP_NAME).bbackz.array`

# Constructor
    VelComp(ynum::Int64, xnum::Int64; znum::Int64=1)

## Arguments
- `ynum::Int64`: Number of grid points in the vertical direction
- `xnum::Int64`: Number of grid points in the horizontal direction
- `znum::Int64=1`: Number of grid points in the z-direction (default: 1 for 2D)
"""
mutable struct VelComp <: AbstractArrayGroup
    btopx::BcArrayFloatState
    btopy::BcArrayFloatState
    btopz::BcArrayFloatState
    bbottomx::BcArrayFloatState
    bbottomy::BcArrayFloatState
    bbottomz::BcArrayFloatState
    bleftx::BcArrayFloatState
    blefty::BcArrayFloatState
    bleftz::BcArrayFloatState
    brightx::BcArrayFloatState
    brighty::BcArrayFloatState
    brightz::BcArrayFloatState
    bfrontx::BcArrayFloatState
    bfronty::BcArrayFloatState
    bfrontz::BcArrayFloatState
    bbackx::BcArrayFloatState
    bbacky::BcArrayFloatState
    bbackz::BcArrayFloatState
end

function VelComp(ynum::Int64, xnum::Int64; znum::Int64=1)::VelComp
    
    return VelComp(
        BcArrayFloatState(
            zeros(Float64, xnum+1, 2), # array
            ADATA.btopx.name, # name
            ADATA.btopx.units, # units
            ADATA.btopx.description
        ),
        BcArrayFloatState(
            zeros(Float64, xnum+1, 2), # array
            ADATA.btopy.name, # name
            ADATA.btopy.units, # units
            ADATA.btopy.description
        ),
        BcArrayFloatState(
            zeros(Float64, xnum+1, 2), # array
            ADATA.btopz.name, # name
            ADATA.btopz.units, # units
            ADATA.btopz.description
        ),
        BcArrayFloatState(
            zeros(Float64, xnum+1, 2), # array
            ADATA.bbottomx.name, # name
            ADATA.bbottomx.units, # units
            ADATA.bbottomx.description
        ),
        BcArrayFloatState(
            zeros(Float64, xnum+1, 2), # array
            ADATA.bbottomy.name, # name
            ADATA.bbottomy.units, # units
            ADATA.bbottomy.description
        ),
        BcArrayFloatState(
            zeros(Float64, xnum+1, 2), # array
            ADATA.bbottomz.name, # name
            ADATA.bbottomz.units, # units
            ADATA.bbottomz.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.bleftx.name, # name
            ADATA.bleftx.units, # units
            ADATA.bleftx.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.blefty.name, # name
            ADATA.blefty.units, # units
            ADATA.blefty.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.bleftz.name, # name
            ADATA.bleftz.units, # units
            ADATA.bleftz.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.brightx.name, # name
            ADATA.brightx.units, # units
            ADATA.brightx.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.brighty.name, # name
            ADATA.brighty.units, # units
            ADATA.brighty.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.brightz.name, # name
            ADATA.brightz.units, # units
            ADATA.brightz.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.bfrontx.name, # name
            ADATA.bfrontx.units, # units
            ADATA.bfrontx.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.bfronty.name, # name
            ADATA.bfronty.units, # units
            ADATA.bfronty.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.bfrontz.name, # name
            ADATA.bfrontz.units, # units
            ADATA.bfrontz.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.bbackx.name, # name
            ADATA.bbackx.units, # units
            ADATA.bbackx.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.bbacky.name, # name
            ADATA.bbacky.units, # units
            ADATA.bbacky.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum+1, 2), # array
            ADATA.bbackz.name, # name
            ADATA.bbackz.units, # units
            ADATA.bbackz.description
        )
    )
end

end # module 