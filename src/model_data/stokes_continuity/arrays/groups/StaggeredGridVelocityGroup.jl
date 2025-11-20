module StaggeredGridVelocityGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "staggered_grid_velocity"

const ADATA = get_eb_arrays()

"""
    StaggeredGridVelocity <: AbstractArrayGroup

Array group containing velocity arrays on staggered grid.

# Fields
- `vy1_old::`[`ScalarArray2DState`](@ref): $(ADATA.vy1_old.description)
- `vy1::`[`ScalarArray2DState`](@ref): $(ADATA.vy1.description)
- `vx1_old::`[`ScalarArray2DState`](@ref): $(ADATA.vx1_old.description)
- `vx1::`[`ScalarArray2DState`](@ref): $(ADATA.vx1.description)

# Nested Dot Access
- `vy1_old = $(ROOT_NAME).$(GRP_NAME).vy1_old.array`
- `vy1 = $(ROOT_NAME).$(GRP_NAME).vy1.array`
- `vx1_old = $(ROOT_NAME).$(GRP_NAME).vx1_old.array`
- `vx1 = $(ROOT_NAME).$(GRP_NAME).vx1.array`

# Constructor
    StaggeredGridVelocity(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `StaggeredGridVelocity`: New StaggeredGridVelocity array group with initialized arrays

"""
mutable struct StaggeredGridVelocity <: AbstractArrayGroup
    vy1_old::ScalarArray2DState
    vy1::ScalarArray2DState
    vx1_old::ScalarArray2DState
    vx1::ScalarArray2DState
end

function StaggeredGridVelocity(ynum::Int, xnum::Int)::StaggeredGridVelocity
    return StaggeredGridVelocity(
        ScalarArray2DState(
            ynum,                        # ynum
            xnum,                        # xnum
            ADATA.vy1_old.name,          # name
            ADATA.vy1_old.units,         # units
            ADATA.vy1_old.grid_type,     # grid_type
            ADATA.vy1_old.description    # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.vy1.name,          # name
            ADATA.vy1.units,         # units
            ADATA.vy1.grid_type,     # grid_type
            ADATA.vy1.description    # description
        ),
        ScalarArray2DState(
            ynum,                        # ynum
            xnum,                        # xnum
            ADATA.vx1_old.name,          # name
            ADATA.vx1_old.units,         # units
            ADATA.vx1_old.grid_type,     # grid_type
            ADATA.vx1_old.description    # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.vx1.name,          # name
            ADATA.vx1.units,         # units
            ADATA.vx1.grid_type,     # grid_type
            ADATA.vx1.description    # description
        )
    )
end

end # module
