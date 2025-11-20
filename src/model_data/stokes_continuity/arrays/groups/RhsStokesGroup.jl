module RhsStokesGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.Arrays.ArrayTypes.RhsStokesArray1D: RhsStokesArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "rhs"

const ADATA = get_eb_arrays()

"""
    RhsStokes <: AbstractArrayGroup

Array group containing right-hand side arrays for Stokes equations.

# Fields
- `RX1::`[`ScalarArray2DState`](@ref): $(ADATA.RX1.description)
- `RY1::`[`ScalarArray2DState`](@ref): $(ADATA.RY1.description)
- `RC1::`[`ScalarArray2DState`](@ref): $(ADATA.RC1.description)
- `RHS::`[`RhsStokesArray1DState`](@ref): $(ADATA.RHS.description)

# Nested Dot Access
- `RX1 = $(ROOT_NAME).$(GRP_NAME).RX1.array`
- `RY1 = $(ROOT_NAME).$(GRP_NAME).RY1.array`
- `RC1 = $(ROOT_NAME).$(GRP_NAME).RC1.array`
- `RHS = $(ROOT_NAME).$(GRP_NAME).RHS.array`

# Constructor
    RhsStokes(ynum::Int, xnum::Int)

Create a new RhsStokes array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `RhsStokes`: New RhsStokes array group with initialized arrays

"""
mutable struct RhsStokes <: AbstractArrayGroup
    RX1::ScalarArray2DState
    RY1::ScalarArray2DState
    RC1::ScalarArray2DState
    RHS::RhsStokesArray1DState
end

function RhsStokes(ynum::Int, xnum::Int)::RhsStokes
    return RhsStokes(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.RX1.name,          # name
            ADATA.RX1.units,         # units
            ADATA.RX1.grid_type,     # grid_type
            ADATA.RX1.description    # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.RY1.name,          # name
            ADATA.RY1.units,         # units
            ADATA.RY1.grid_type,     # grid_type
            ADATA.RY1.description    # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.RC1.name,          # name
            ADATA.RC1.units,         # units
            ADATA.RC1.grid_type,     # grid_type
            ADATA.RC1.description    # description
        ),
        RhsStokesArray1DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.RHS.name,          # name
            ADATA.RHS.units,         # units
            ADATA.RHS.description    # description
        )
    )
end

end # module
