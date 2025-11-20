module ResidualStokesGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.Arrays.ArrayTypes.SolutionArray1D: SolutionArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "residuals"

const ADATA = get_eb_arrays()

"""
    ResidualStokes <: AbstractArrayGroup

Array group containing residual arrays for Stokes equations.

# Fields
- `resx1::`[`ScalarArray2DState`](@ref): $(ADATA.resx1.description)
- `resnlx::`[`ScalarArray2DState`](@ref): $(ADATA.resnlx.description)
- `resy1::`[`ScalarArray2DState`](@ref): $(ADATA.resy1.description)
- `resnly::`[`ScalarArray2DState`](@ref): $(ADATA.resnly.description)
- `resc1::`[`ScalarArray2DState`](@ref): $(ADATA.resc1.description)
- `resnlc::`[`ScalarArray2DState`](@ref): $(ADATA.resnlc.description)
- `resnl::`[`SolutionArray1DState`](@ref): $(ADATA.resnl.description)

# Nested Dot Access
- `resx1 = $(ROOT_NAME).$(GRP_NAME).resx1.array`
- `resy1 = $(ROOT_NAME).$(GRP_NAME).resy1.array`
- `resc1 = $(ROOT_NAME).$(GRP_NAME).resc1.array`
- `resnl = $(ROOT_NAME).$(GRP_NAME).resnl.array`

# Constructor
    ResidualStokes(ynum::Int, xnum::Int)

Create a new ResidualStokes array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `ResidualStokes`: New ResidualStokes array group with initialized arrays

"""
mutable struct ResidualStokes <: AbstractArrayGroup
    resx1::ScalarArray2DState
    resnlx::ScalarArray2DState
    resy1::ScalarArray2DState
    resnly::ScalarArray2DState
    resc1::ScalarArray2DState
    resnlc::ScalarArray2DState
    resnl::SolutionArray1DState
end

function ResidualStokes(ynum::Int, xnum::Int)::ResidualStokes
    return ResidualStokes(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.resx1.name,        # name
            ADATA.resx1.units,       # units
            ADATA.resx1.grid_type,   # grid_type
            ADATA.resx1.description  # description
        ),
        ScalarArray2DState(
            ynum,                     # ynum
            xnum,                     # xnum
            ADATA.resnlx.name,        # name
            ADATA.resnlx.units,       # units
            ADATA.resnlx.grid_type,   # grid_type
            ADATA.resnlx.description  # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.resy1.name,        # name
            ADATA.resy1.units,       # units
            ADATA.resy1.grid_type,   # grid_type
            ADATA.resy1.description  # description
        ),
        ScalarArray2DState(
            ynum,                     # ynum
            xnum,                     # xnum
            ADATA.resnly.name,        # name
            ADATA.resnly.units,       # units
            ADATA.resnly.grid_type,   # grid_type
            ADATA.resnly.description  # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.resc1.name,        # name
            ADATA.resc1.units,       # units
            ADATA.resc1.grid_type,   # grid_type
            ADATA.resc1.description  # description
        ),
        ScalarArray2DState(
            ynum,                     # ynum
            xnum,                     # xnum
            ADATA.resnlc.name,        # name
            ADATA.resnlc.units,       # units
            ADATA.resnlc.grid_type,   # grid_type
            ADATA.resnlc.description  # description
        ),
        SolutionArray1DState(
            ynum,                        # ynum
            xnum,                        # xnum
            ADATA.resnl.name,            # name
            ADATA.resnl.units,           # units
            ADATA.resnl.description,     # description
            "normal"                     # array_type
        )
    )
end

end # module
