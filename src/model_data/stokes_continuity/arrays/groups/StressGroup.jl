module StressGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "stress"

const ADATA = get_eb_arrays()

"""
    Stress <: AbstractArrayGroup

Array group containing stress tensor component arrays.

# Fields
- `sxy0::`[`ScalarArray2DState`](@ref): $(ADATA.sxy0.description)
- `sxy1::`[`ScalarArray2DState`](@ref): $(ADATA.sxy1.description)
- `sxy2::`[`ScalarArray2DState`](@ref): $(ADATA.sxy2.description)
- `sxx0::`[`ScalarArray2DState`](@ref): $(ADATA.sxx0.description)
- `sxx1::`[`ScalarArray2DState`](@ref): $(ADATA.sxx1.description)
- `sxx2::`[`ScalarArray2DState`](@ref): $(ADATA.sxx2.description)

# Nested Dot Access
- `sxy0 = $(ROOT_NAME).$(GRP_NAME).sxy0.array`
- `sxy1 = $(ROOT_NAME).$(GRP_NAME).sxy1.array`
- `sxy2 = $(ROOT_NAME).$(GRP_NAME).sxy2.array`
- `sxx0 = $(ROOT_NAME).$(GRP_NAME).sxx0.array`

# Constructor
    Stress(ynum::Int, xnum::Int)

Create a new Stress array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `Stress`: New Stress array group with initialized arrays

"""
mutable struct Stress <: AbstractArrayGroup
    sxy0::ScalarArray2DState
    sxy1::ScalarArray2DState
    sxy2::ScalarArray2DState
    sxx0::ScalarArray2DState
    sxx1::ScalarArray2DState
    sxx2::ScalarArray2DState
end

function Stress(ynum::Int, xnum::Int)::Stress
    data = Stress(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.sxy0.name,         # name
            ADATA.sxy0.units,        # units
            ADATA.sxy0.grid_type,    # grid_type
            ADATA.sxy0.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.sxy1.name,         # name
            ADATA.sxy1.units,        # units
            ADATA.sxy1.grid_type,    # grid_type
            ADATA.sxy1.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.sxy2.name,         # name
            ADATA.sxy2.units,        # units
            ADATA.sxy2.grid_type,    # grid_type
            ADATA.sxy2.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.sxx0.name,         # name
            ADATA.sxx0.units,        # units
            ADATA.sxx0.grid_type,    # grid_type
            ADATA.sxx0.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.sxx1.name,         # name
            ADATA.sxx1.units,        # units
            ADATA.sxx1.grid_type,    # grid_type
            ADATA.sxx1.description   # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.sxx2.name,         # name
            ADATA.sxx2.units,        # units
            ADATA.sxx2.grid_type,    # grid_type
            ADATA.sxx2.description   # description
        )
    )
    update_output_format(data)
    return data
end

function update_output_format(data::Stress)::Nothing
    data.sxy0.outform.fac1 = 1e-6
    data.sxy0.outform.units = "MPa"
    data.sxy0.outform.fname = "Sxy0_MPa"
    data.sxy1.outform.fac1 = 1e-6
    data.sxy1.outform.units = "MPa"
    data.sxy1.outform.fname = "Sxy1_MPa"
    data.sxy2.outform.fac1 = 1e-6
    data.sxy2.outform.units = "MPa"
    data.sxy2.outform.fname = "Sxy_MPa"
    data.sxx0.outform.fac1 = 1e-6
    data.sxx0.outform.units = "MPa"
    data.sxx0.outform.fname = "Sxx0_MPa"
    data.sxx1.outform.fac1 = 1e-6
    data.sxx1.outform.units = "MPa"
    data.sxx1.outform.fname = "Sxx1_MPa"
    data.sxx2.outform.fac1 = 1e-6
    data.sxx2.outform.units = "MPa"
    data.sxx2.outform.fname = "Sxx_MPa"
    return nothing
end

end # module
