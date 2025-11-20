module ShearModulusGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "shear_modulus"

const ADATA = get_eb_arrays()

"""
    ShearModulus <: AbstractArrayGroup

Array group containing shear modulus arrays.

# Fields
- `mus0::`[`ScalarArray2DState`](@ref): $(ADATA.mus0.description)
- `mus1::`[`ScalarArray2DState`](@ref): $(ADATA.mus1.description)
- `mun0::`[`ScalarArray2DState`](@ref): $(ADATA.mun0.description)
- `mun1::`[`ScalarArray2DState`](@ref): $(ADATA.mun1.description)

# Nested Dot Access
- `mus0 = $(ROOT_NAME).$(GRP_NAME).mus0.array`
- `mus1 = $(ROOT_NAME).$(GRP_NAME).mus1.array`
- `mun0 = $(ROOT_NAME).$(GRP_NAME).mun0.array`
- `mun1 = $(ROOT_NAME).$(GRP_NAME).mun1.array`

# Constructor
    ShearModulus(ynum::Int, xnum::Int)

Create a new ShearModulus array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `ShearModulus`: New ShearModulus array group with initialized arrays

"""
mutable struct ShearModulus <: AbstractArrayGroup
    mus0::ScalarArray2DState
    mus1::ScalarArray2DState
    mun0::ScalarArray2DState
    mun1::ScalarArray2DState
end

function ShearModulus(ynum::Int, xnum::Int)::ShearModulus
    return ShearModulus(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.mus0.name,         # name
            ADATA.mus0.units,        # units
            ADATA.mus0.grid_type,    # grid_type
            ADATA.mus0.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.mus1.name,         # name
            ADATA.mus1.units,        # units
            ADATA.mus1.grid_type,    # grid_type
            ADATA.mus1.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.mun0.name,         # name
            ADATA.mun0.units,        # units
            ADATA.mun0.grid_type,    # grid_type
            ADATA.mun0.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.mun1.name,         # name
            ADATA.mun1.units,        # units
            ADATA.mun1.grid_type,    # grid_type
            ADATA.mun1.description   # description
        )
    )
end

end # module
