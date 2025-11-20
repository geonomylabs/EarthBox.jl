module RhsHeatGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.Arrays.ArrayTypes.RhsHeatArray1D: RhsHeatArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.heat_equation.arrays"
const GRP_NAME = "rhs"

const ADATA = get_eb_arrays()

"""
    RhsHeat <: AbstractArrayGroup

Array group containing right-hand side arrays for heat equation.

# Fields
- `RT1::`[`ScalarArray2DState`](@ref): $(ADATA.RT1.description)
- `RHSheat::`[`RhsHeatArray1DState`](@ref): $(ADATA.RHSheat.description)

# Nested Dot Access
- `RT1 = $(ROOT_NAME).$(GRP_NAME).RT1.array`
- `RHSheat = $(ROOT_NAME).$(GRP_NAME).RHSheat.array`

# Constructor
    RhsHeat(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `RhsHeat`: New RhsHeat array group with initialized arrays

"""
mutable struct RhsHeat <: AbstractArrayGroup
    RT1::ScalarArray2DState
    RHSheat::RhsHeatArray1DState
end

function RhsHeat(ynum::Int, xnum::Int)::RhsHeat
    return RhsHeat(
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.RT1.name,              # name
            ADATA.RT1.units,             # units
            ADATA.RT1.grid_type,         # grid_type
            ADATA.RT1.description        # description
        ),
        RhsHeatArray1DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.RHSheat.name,          # name
            ADATA.RHSheat.units,         # units
            ADATA.RHSheat.description    # description
        )
    )
end

end # module
