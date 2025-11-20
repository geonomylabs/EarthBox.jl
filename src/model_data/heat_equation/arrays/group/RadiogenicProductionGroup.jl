module RadiogenicProductionGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.heat_equation.arrays"
const GRP_NAME = "radiogenic_production"

const ADATA = get_eb_arrays()

"""
    RadiogenicProduction <: AbstractArrayGroup

Array group containing radiogenic heat production arrays.

# Fields
- `hr0::`[`ScalarArray2DState`](@ref): $(ADATA.hr0.description)
- `hr1::`[`ScalarArray2DState`](@ref): $(ADATA.hr1.description)

# Nested Dot Access
- `hr0 = $(ROOT_NAME).$(GRP_NAME).hr0.array`
- `hr1 = $(ROOT_NAME).$(GRP_NAME).hr1.array`

# Constructor
    RadiogenicProduction(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `RadiogenicProduction`: New RadiogenicProduction array group with initialized arrays

"""
mutable struct RadiogenicProduction <: AbstractArrayGroup
    hr0::ScalarArray2DState
    hr1::ScalarArray2DState
end

function RadiogenicProduction(ynum::Int, xnum::Int)::RadiogenicProduction
    return RadiogenicProduction(
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.hr0.name,              # name
            ADATA.hr0.units,             # units
            ADATA.hr0.grid_type,         # grid_type
            ADATA.hr0.description        # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.hr1.name,              # name
            ADATA.hr1.units,             # units
            ADATA.hr1.grid_type,         # grid_type
            ADATA.hr1.description        # description
        )
    )
end

end # module
