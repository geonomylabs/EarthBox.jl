module ThermalConductivityGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.heat_equation.arrays"
const GRP_NAME = "thermal_conductivity"

const ADATA = get_eb_arrays()

"""
    ThermalConductivity <: AbstractArrayGroup

Array group containing thermal conductivity arrays.

# Fields
- `kt0::`[`ScalarArray2DState`](@ref): $(ADATA.kt0.description)
- `kt1::`[`ScalarArray2DState`](@ref): $(ADATA.kt1.description)

# Nested Dot Access
- `kt0 = $(ROOT_NAME).$(GRP_NAME).kt0.array`
- `kt1 = $(ROOT_NAME).$(GRP_NAME).kt1.array`

# Constructor
    ThermalConductivity(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `ThermalConductivity`: New ThermalConductivity array group with initialized arrays

"""
mutable struct ThermalConductivity <: AbstractArrayGroup
    kt0::ScalarArray2DState
    kt1::ScalarArray2DState
end

function ThermalConductivity(ynum::Int, xnum::Int)::ThermalConductivity
    state = ThermalConductivity(
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.kt0.name,              # name
            ADATA.kt0.units,             # units
            ADATA.kt0.grid_type,         # grid_type
            ADATA.kt0.description        # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.kt1.name,              # name
            ADATA.kt1.units,             # units
            ADATA.kt1.grid_type,         # grid_type
            ADATA.kt1.description        # description
        )
    )
    update_output_format!(state)
    return state
end

function update_output_format!(state::ThermalConductivity)::Nothing
    state.kt0.outform.fname = "therm_cond0"
    state.kt1.outform.fname = "therm_cond"
    return nothing
end

end # module
