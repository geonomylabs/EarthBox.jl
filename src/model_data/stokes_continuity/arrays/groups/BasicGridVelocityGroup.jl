module BasicGridVelocityGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "basic_grid_velocity"

const ADATA = get_eb_arrays()

"""
    BasicGridVelocity <: AbstractArrayGroup

Array group containing velocity arrays on basic grid.

# Fields
- `vxb::`[`ScalarArray2DState`](@ref): $(ADATA.vxb.description)
- `vyb::`[`ScalarArray2DState`](@ref): $(ADATA.vyb.description)

# Nested Dot Access
- `vxb = $(ROOT_NAME).$(GRP_NAME).vxb.array`
- `vyb = $(ROOT_NAME).$(GRP_NAME).vyb.array`

# Constructor
    BasicGridVelocity(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `BasicGridVelocity`: New BasicGridVelocity array group with initialized arrays

"""
mutable struct BasicGridVelocity <: AbstractArrayGroup
    vxb::ScalarArray2DState
    vyb::ScalarArray2DState
end

function BasicGridVelocity(ynum::Int, xnum::Int)::BasicGridVelocity
    state = BasicGridVelocity(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.vxb.name,          # name
            ADATA.vxb.units,         # units
            ADATA.vxb.grid_type,     # grid_type
            ADATA.vxb.description    # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.vyb.name,          # name
            ADATA.vyb.units,         # units
            ADATA.vyb.grid_type,     # grid_type
            ADATA.vyb.description    # description
        )
    )
    update_output_format(state)
    return state
end

function update_output_format(data::BasicGridVelocity)::Nothing
    conversion_factor = 100.0 * 365.25 * 24.0 * 3600.0
    
    data.vxb.outform.fac1 = conversion_factor
    data.vxb.outform.units = "cm/yr"
    data.vxb.outform.fname = "vx"
    data.vxb.outform.header = "VX_(cm/yr)"
    
    data.vyb.outform.fac1 = conversion_factor
    data.vyb.outform.units = "cm/yr"
    data.vyb.outform.fname = "vy"
    data.vyb.outform.header = "VY_(cm/yr)"
    
    return nothing
end

end # module
