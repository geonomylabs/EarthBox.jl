module AdiabaticProductionGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.heat_equation.arrays"
const GRP_NAME = "adiabatic_production"

const ADATA = get_eb_arrays()

"""
    AdiabaticProduction <: AbstractArrayGroup

Array group containing adiabatic heat production arrays.

# Fields
- `ha0::`[`ScalarArray2DState`](@ref): $(ADATA.ha0.description)
- `ha1::`[`ScalarArray2DState`](@ref): $(ADATA.ha1.description)

# Nested Dot Access
- `ha0 = $(ROOT_NAME).$(GRP_NAME).ha0.array`
- `ha1 = $(ROOT_NAME).$(GRP_NAME).ha1.array`

# Constructor
    AdiabaticProduction(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `AdiabaticProduction`: New AdiabaticProduction array group with initialized arrays

"""
mutable struct AdiabaticProduction <: AbstractArrayGroup
    ha0::ScalarArray2DState
    ha1::ScalarArray2DState
end

function AdiabaticProduction(ynum::Int, xnum::Int)::AdiabaticProduction
    return AdiabaticProduction(
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.ha0.name,              # name
            ADATA.ha0.units,             # units
            ADATA.ha0.grid_type,         # grid_type
            ADATA.ha0.description        # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.ha1.name,              # name
            ADATA.ha1.units,             # units
            ADATA.ha1.grid_type,         # grid_type
            ADATA.ha1.description        # description
        )
    )
end

end # module
