module DensityGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "density"

const ADATA = get_eb_arrays()

"""
    Density <: AbstractArrayGroup

Array group containing density field arrays.

# Fields
- `rho0::`[`ScalarArray2DState`](@ref): $(ADATA.rho0.description)
- `rho1::`[`ScalarArray2DState`](@ref): $(ADATA.rho1.description)
- `rho0_vy::`[`ScalarArray2DState`](@ref): $(ADATA.rho0_vy.description)
- `rho1_vy::`[`ScalarArray2DState`](@ref): $(ADATA.rho1_vy.description)

# Nested Dot Access
- `rho0 = $(ROOT_NAME).$(GRP_NAME).rho0.array`
- `rho1 = $(ROOT_NAME).$(GRP_NAME).rho1.array`
- `rho0_vy = $(ROOT_NAME).$(GRP_NAME).rho0_vy.array`
- `rho1_vy = $(ROOT_NAME).$(GRP_NAME).rho1_vy.array`

# Constructor
    Density(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `Density`: New Density array group with initialized arrays

"""
mutable struct Density <: AbstractArrayGroup
    rho0::ScalarArray2DState
    rho1::ScalarArray2DState
    rho0_vy::ScalarArray2DState
    rho1_vy::ScalarArray2DState
end

function Density(ynum::Int, xnum::Int)::Density
    state = Density(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.rho0.name,         # name
            ADATA.rho0.units,        # units
            ADATA.rho0.grid_type,    # grid_type
            ADATA.rho0.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.rho1.name,         # name
            ADATA.rho1.units,        # units
            ADATA.rho1.grid_type,    # grid_type
            ADATA.rho1.description   # description
        ),
        ScalarArray2DState(
            ynum,                     # ynum
            xnum,                     # xnum
            ADATA.rho0_vy.name,       # name
            ADATA.rho0_vy.units,      # units
            ADATA.rho0_vy.grid_type,  # grid_type
            ADATA.rho0_vy.description # description
        ),
        ScalarArray2DState(
            ynum,                     # ynum
            xnum,                     # xnum
            ADATA.rho1_vy.name,       # name
            ADATA.rho1_vy.units,      # units
            ADATA.rho1_vy.grid_type,  # grid_type
            ADATA.rho1_vy.description # description
        )
    )
    update_output_format(state)
    return state
end

function update_output_format(data::Density)::Nothing
    data.rho0.outform.fname = "rho0_kg_m3"
    data.rho1.outform.fname = "rho_kg_m3"
    data.rho0_vy.outform.fname = "rho0_vy_kg_m3"
    data.rho1_vy.outform.fname = "rho_vy_kg_m3"
    return nothing
end

end # module
