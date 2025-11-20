"""
Module for temperature boundary condition arrays.

Provides data structures for specifying temperature boundary conditions along
the domain boundaries.
"""
module TemperatureGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.EarthBoxDtypes: AbstractArrayGroup
import EarthBox.Arrays.ArrayTypes.BcArrayFloat: BcArrayFloatState

const ROOT_NAME = "model.bcs.arrays"
const GRP_NAME = "temperature"

const ADATA = get_eb_arrays()

"""
    Temperature <: AbstractArrayGroup

Array group for temperature boundary conditions.

# Fields
- `btopt::`[`BcArrayFloatState`](@ref): $(ADATA.btopt.description)
- `bbottomt::`[`BcArrayFloatState`](@ref): $(ADATA.bbottomt.description)
- `bleftt::`[`BcArrayFloatState`](@ref): $(ADATA.bleftt.description)
- `brightt::`[`BcArrayFloatState`](@ref): $(ADATA.brightt.description)

# Constructor
    Temperature(ynum::Int, xnum::Int)

## Arguments
- `ynum::Int`: Number of grid points in the vertical direction
- `xnum::Int`: Number of grid points in the horizontal direction
"""
mutable struct Temperature <: AbstractArrayGroup
    btopt::BcArrayFloatState
    bbottomt::BcArrayFloatState
    bleftt::BcArrayFloatState
    brightt::BcArrayFloatState
end

function Temperature(ynum::Int, xnum::Int)::Temperature
    
    return Temperature(
        BcArrayFloatState(
            zeros(Float64, xnum, 2), # array
            ADATA.btopt.name, # name
            ADATA.btopt.units, # units
            ADATA.btopt.description
        ),
        BcArrayFloatState(
            zeros(Float64, xnum, 2), # array
            ADATA.bbottomt.name, # name
            ADATA.bbottomt.units, # units
            ADATA.bbottomt.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum, 2), # array
            ADATA.bleftt.name, # name
            ADATA.bleftt.units, # units
            ADATA.bleftt.description
        ),
        BcArrayFloatState(
            zeros(Float64, ynum, 2), # array
            ADATA.brightt.name, # name
            ADATA.brightt.units, # units
            ADATA.brightt.description
        )
    )
end

end # module 