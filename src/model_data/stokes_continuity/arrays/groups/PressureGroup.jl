module PressureGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "pressure"

const ADATA = get_eb_arrays()

"""
    Pressure <: AbstractArrayGroup

Array group containing pressure field arrays.

# Fields
- `pr1_old::`[`ScalarArray2DState`](@ref): $(ADATA.pr1_old.description)
- `pr1::`[`ScalarArray2DState`](@ref): $(ADATA.pr1.description)

# Nested Dot Access
- `pr1_old = $(ROOT_NAME).$(GRP_NAME).pr1_old.array`
- `pr1 = $(ROOT_NAME).$(GRP_NAME).pr1.array`

# Constructor
    Pressure(ynum::Int, xnum::Int)

Create a new Pressure array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `Pressure`: New Pressure array group with initialized arrays

"""
mutable struct Pressure <: AbstractArrayGroup
    pr1_old::ScalarArray2DState
    pr1::ScalarArray2DState
end

function Pressure(ynum::Int, xnum::Int)::Pressure
    state = Pressure(
        ScalarArray2DState(
            ynum,                      # ynum
            xnum,                      # xnum
            ADATA.pr1_old.name,        # name
            ADATA.pr1_old.units,       # units
            ADATA.pr1_old.grid_type,   # grid_type
            ADATA.pr1_old.description  # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.pr1.name,          # name
            ADATA.pr1.units,         # units
            ADATA.pr1.grid_type,     # grid_type
            ADATA.pr1.description    # description
        )
    )
    update_output_format(state)
    return state
end

function update_output_format(data::Pressure)::Nothing
    data.pr1_old.outform.fac1 = 1e-9
    data.pr1_old.outform.units = "GPa"
    data.pr1_old.outform.fname = "pressure_old_GPa"
    
    data.pr1.outform.fac1 = 1e-9
    data.pr1.outform.units = "GPa"
    data.pr1.outform.fname = "pressure_GPa"
    return nothing
end

end # module
