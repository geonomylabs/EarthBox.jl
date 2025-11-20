module TemperatureGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.heat_equation.arrays"
const GRP_NAME = "temperature"

const ADATA = get_eb_arrays()

"""
    Temperature <: AbstractArrayGroup

Array group containing temperature field arrays.

# Fields
- `tk0::`[`ScalarArray2DState`](@ref): $(ADATA.tk0.description)
- `tk1::`[`ScalarArray2DState`](@ref): $(ADATA.tk1.description)
- `tk2::`[`ScalarArray2DState`](@ref): $(ADATA.tk2.description)
- `dtk1::`[`ScalarArray2DState`](@ref): $(ADATA.dtk1.description)
- `dtkn::`[`ScalarArray2DState`](@ref): $(ADATA.dtkn.description)

# Nested Dot Access
- `tk0 = $(ROOT_NAME).$(GRP_NAME).tk0.array`
- `tk1 = $(ROOT_NAME).$(GRP_NAME).tk1.array`
- `tk2 = $(ROOT_NAME).$(GRP_NAME).tk2.array`
- `dtk1 = $(ROOT_NAME).$(GRP_NAME).dtk1.array`
- `dtkn = $(ROOT_NAME).$(GRP_NAME).dtkn.array`

# Constructor
    Temperature(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `Temperature`: New Temperature array group with initialized arrays

"""
mutable struct Temperature <: AbstractArrayGroup
    tk0::ScalarArray2DState
    tk1::ScalarArray2DState
    tk2::ScalarArray2DState
    dtk1::ScalarArray2DState
    dtkn::ScalarArray2DState
end

function Temperature(ynum::Int, xnum::Int)::Temperature
    state = Temperature(
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.tk0.name,              # name
            ADATA.tk0.units,             # units
            ADATA.tk0.grid_type,         # grid_type
            ADATA.tk0.description        # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.tk1.name,              # name
            ADATA.tk1.units,             # units
            ADATA.tk1.grid_type,         # grid_type
            ADATA.tk1.description        # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.tk2.name,              # name
            ADATA.tk2.units,             # units
            ADATA.tk2.grid_type,         # grid_type
            ADATA.tk2.description        # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.dtk1.name,             # name
            ADATA.dtk1.units,            # units
            ADATA.dtk1.grid_type,        # grid_type
            ADATA.dtk1.description       # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.dtkn.name,             # name
            ADATA.dtkn.units,            # units
            ADATA.dtkn.grid_type,        # grid_type
            ADATA.dtkn.description       # description
        )
    )
    update_output_format!(state)
    return state
end

function update_output_format!(state::Temperature)::Nothing
    state.tk0.outform.fac2 = -273.0
    state.tk0.outform.units = "C"
    state.tk0.outform.fname = "OldTempC"
    
    state.tk1.outform.fac2 = -273.0
    state.tk1.outform.units = "C"
    state.tk1.outform.fname = "TempC"
    
    state.tk2.outform.fac2 = -273.0
    state.tk2.outform.units = "C"
    state.tk2.outform.fname = "ConductiveTempC"
    
    state.dtk1.outform.fac2 = -273.0
    state.dtk1.outform.units = "C"
    state.dtk1.outform.fname = "DeltaGridTempC"
    
    state.dtkn.outform.fac2 = -273.0
    state.dtkn.outform.units = "C"
    state.dtkn.outform.fname = "SubgridDeltaTempC"
    
    return nothing
end

end # module
