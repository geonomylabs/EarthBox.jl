module HeatResidualGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.heat_equation.arrays"
const GRP_NAME = "residuals"

const ADATA = get_eb_arrays()

"""
    HeatResidual <: AbstractArrayGroup

Array group containing heat equation residual array.

# Fields
- `rest::`[`ScalarArray2DState`](@ref): $(ADATA.rest.description)

# Nested Dot Access
- `rest = $(ROOT_NAME).$(GRP_NAME).rest.array`

# Constructor
    HeatResidual(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `HeatResidual`: New HeatResidual array group with initialized arrays

"""
mutable struct HeatResidual <: AbstractArrayGroup
    rest::ScalarArray2DState
end

function HeatResidual(ynum::Int, xnum::Int)::HeatResidual
    return HeatResidual(
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.rest.name,              # name
            ADATA.rest.units,             # units
            ADATA.rest.grid_type,         # grid_type
            ADATA.rest.description        # description
        )
    )
end

end # module
