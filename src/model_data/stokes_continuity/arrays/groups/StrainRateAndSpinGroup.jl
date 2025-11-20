module StrainRateAndSpinGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "strain_rate_and_spin"

const ADATA = get_eb_arrays()

"""
    StrainRateAndSpin <: AbstractArrayGroup

Array group containing strain rate and spin tensor arrays.

# Fields
- `exy::`[`ScalarArray2DState`](@ref): $(ADATA.exy.description)
- `exx::`[`ScalarArray2DState`](@ref): $(ADATA.exx.description)
- `eii::`[`ScalarArray2DState`](@ref): $(ADATA.eii.description)
- `esp::`[`ScalarArray2DState`](@ref): $(ADATA.esp.description)
- `eii_plastic_basic::`[`ScalarArray2DState`](@ref): $(ADATA.eii_plastic_basic.description)
- `eii_plastic_pressure::`[`ScalarArray2DState`](@ref): $(ADATA.eii_plastic_pressure.description)

# Nested Dot Access
- `exy = $(ROOT_NAME).$(GRP_NAME).exy.array`
- `exx = $(ROOT_NAME).$(GRP_NAME).exx.array`
- `eii = $(ROOT_NAME).$(GRP_NAME).eii.array`
- `esp = $(ROOT_NAME).$(GRP_NAME).esp.array`

# Constructor
    StrainRateAndSpin(ynum::Int, xnum::Int)

Create a new StrainRateAndSpin array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `StrainRateAndSpin`: New StrainRateAndSpin array group with initialized arrays

"""
mutable struct StrainRateAndSpin <: AbstractArrayGroup
    exy::ScalarArray2DState
    exx::ScalarArray2DState
    eii::ScalarArray2DState
    esp::ScalarArray2DState
    eii_plastic_basic::ScalarArray2DState
    eii_plastic_pressure::ScalarArray2DState
end

function StrainRateAndSpin(ynum::Int, xnum::Int)::StrainRateAndSpin
    data = StrainRateAndSpin(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.exy.name,          # name
            ADATA.exy.units,         # units
            ADATA.exy.grid_type,     # grid_type
            ADATA.exy.description    # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.exx.name,          # name
            ADATA.exx.units,         # units
            ADATA.exx.grid_type,     # grid_type
            ADATA.exx.description    # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.eii.name,          # name
            ADATA.eii.units,         # units
            ADATA.eii.grid_type,     # grid_type
            ADATA.eii.description    # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.esp.name,          # name
            ADATA.esp.units,         # units
            ADATA.esp.grid_type,     # grid_type
            ADATA.esp.description    # description
        ),
        ScalarArray2DState(
            ynum,                              # ynum
            xnum,                              # xnum
            ADATA.eii_plastic_basic.name,      # name
            ADATA.eii_plastic_basic.units,     # units
            ADATA.eii_plastic_basic.grid_type, # grid_type
            ADATA.eii_plastic_basic.description # description
        ),
        ScalarArray2DState(
            ynum,                                 # ynum
            xnum,                                 # xnum
            ADATA.eii_plastic_pressure.name,      # name
            ADATA.eii_plastic_pressure.units,     # units
            ADATA.eii_plastic_pressure.grid_type, # grid_type
            ADATA.eii_plastic_pressure.description # description
        )
    )
    update_output_format(data)
    return data
end

function update_output_format(data::StrainRateAndSpin)::Nothing
    data.exy.outform.log10 = true
    data.exy.outform.fname = "Exy_log"
    data.exx.outform.log10 = true
    data.exx.outform.fname = "Exx_log"
    data.eii.outform.log10 = true
    data.eii.outform.fname = "Eii_log"
    data.esp.outform.log10 = true
    data.esp.outform.fname = "Esp_log"
    return nothing
end

end # module
