module StressChangeGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "stress_change"

const ADATA = get_eb_arrays()

"""
    StressChange <: AbstractArrayGroup

Array group containing stress change arrays.

# Fields
- `dsxy::`[`ScalarArray2DState`](@ref): $(ADATA.dsxy.description)
- `dsxx::`[`ScalarArray2DState`](@ref): $(ADATA.dsxx.description)
- `dsxyn::`[`ScalarArray2DState`](@ref): $(ADATA.dsxyn.description)
- `dsxxn::`[`ScalarArray2DState`](@ref): $(ADATA.dsxxn.description)

# Nested Dot Access
- `dsxy = $(ROOT_NAME).$(GRP_NAME).dsxy.array`
- `dsxx = $(ROOT_NAME).$(GRP_NAME).dsxx.array`
- `dsxyn = $(ROOT_NAME).$(GRP_NAME).dsxyn.array`
- `dsxxn = $(ROOT_NAME).$(GRP_NAME).dsxxn.array`

# Constructor
    StressChange(ynum::Int, xnum::Int)

Create a new StressChange array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `StressChange`: New StressChange array group with initialized arrays

"""
mutable struct StressChange <: AbstractArrayGroup
    dsxy::ScalarArray2DState
    dsxx::ScalarArray2DState
    dsxyn::ScalarArray2DState
    dsxxn::ScalarArray2DState
end

function StressChange(ynum::Int, xnum::Int)::StressChange
    return StressChange(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.dsxy.name,         # name
            ADATA.dsxy.units,        # units
            ADATA.dsxy.grid_type,    # grid_type
            ADATA.dsxy.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.dsxx.name,         # name
            ADATA.dsxx.units,        # units
            ADATA.dsxx.grid_type,    # grid_type
            ADATA.dsxx.description   # description
        ),
        ScalarArray2DState(
            ynum,                     # ynum
            xnum,                     # xnum
            ADATA.dsxyn.name,         # name
            ADATA.dsxyn.units,        # units
            ADATA.dsxyn.grid_type,    # grid_type
            ADATA.dsxyn.description   # description
        ),
        ScalarArray2DState(
            ynum,                     # ynum
            xnum,                     # xnum
            ADATA.dsxxn.name,         # name
            ADATA.dsxxn.units,        # units
            ADATA.dsxxn.grid_type,    # grid_type
            ADATA.dsxxn.description   # description
        )
    )
end

end # module
