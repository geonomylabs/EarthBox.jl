module RhoCpGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.heat_equation.arrays"
const GRP_NAME = "rhocp"

const ADATA = get_eb_arrays()

"""
    RhoCp <: AbstractArrayGroup

Array group containing density times heat capacity arrays.

# Fields
- `rhocp0::`[`ScalarArray2DState`](@ref): $(ADATA.rhocp0.description)
- `rhocp1::`[`ScalarArray2DState`](@ref): $(ADATA.rhocp1.description)

# Nested Dot Access
- `rhocp0 = $(ROOT_NAME).$(GRP_NAME).rhocp0.array`
- `rhocp1 = $(ROOT_NAME).$(GRP_NAME).rhocp1.array`

# Constructor
    RhoCp(ynum::Int, xnum::Int)

Create a new RhoCp array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `RhoCp`: New RhoCp array group with initialized arrays

"""
mutable struct RhoCp <: AbstractArrayGroup
    rhocp0::ScalarArray2DState
    rhocp1::ScalarArray2DState
end

function RhoCp(ynum::Int, xnum::Int)::RhoCp
    return RhoCp(
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.rhocp0.name,           # name
            ADATA.rhocp0.units,          # units
            ADATA.rhocp0.grid_type,      # grid_type
            ADATA.rhocp0.description     # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.rhocp1.name,           # name
            ADATA.rhocp1.units,          # units
            ADATA.rhocp1.grid_type,      # grid_type
            ADATA.rhocp1.description     # description
        )
    )
end

end # module
