module ViscosityGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "viscosity"

const ADATA = get_eb_arrays()

"""
    Viscosity <: AbstractArrayGroup

Array group containing viscosity field arrays.

# Fields
- `etas0::`[`ScalarArray2DState`](@ref): $(ADATA.etas0.description)
- `etas1::`[`ScalarArray2DState`](@ref): $(ADATA.etas1.description)
- `etan0::`[`ScalarArray2DState`](@ref): $(ADATA.etan0.description)
- `etan1::`[`ScalarArray2DState`](@ref): $(ADATA.etan1.description)
- `eta_flow::`[`ScalarArray2DState`](@ref): $(ADATA.eta_flow.description)
- `eta_flow0::`[`ScalarArray2DState`](@ref): $(ADATA.eta_flow0.description)

# Nested Dot Access
- `etas0 = $(ROOT_NAME).$(GRP_NAME).etas0.array`
- `etas1 = $(ROOT_NAME).$(GRP_NAME).etas1.array`
- `etan0 = $(ROOT_NAME).$(GRP_NAME).etan0.array`
- `etan1 = $(ROOT_NAME).$(GRP_NAME).etan1.array`

# Constructor
    Viscosity(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `Viscosity`: New Viscosity array group with initialized arrays

"""
mutable struct Viscosity <: AbstractArrayGroup
    etas0::ScalarArray2DState
    etas1::ScalarArray2DState
    etan0::ScalarArray2DState
    etan1::ScalarArray2DState
    eta_flow::ScalarArray2DState
    eta_flow0::ScalarArray2DState
end

function Viscosity(ynum::Int, xnum::Int)::Viscosity
    data = Viscosity(
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.etas0.name,        # name
            ADATA.etas0.units,       # units
            ADATA.etas0.grid_type,   # grid_type
            ADATA.etas0.description  # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.etas1.name,        # name
            ADATA.etas1.units,       # units
            ADATA.etas1.grid_type,   # grid_type
            ADATA.etas1.description  # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.etan0.name,        # name
            ADATA.etan0.units,       # units
            ADATA.etan0.grid_type,   # grid_type
            ADATA.etan0.description  # description
        ),
        ScalarArray2DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.etan1.name,        # name
            ADATA.etan1.units,       # units
            ADATA.etan1.grid_type,   # grid_type
            ADATA.etan1.description  # description
        ),
        ScalarArray2DState(
            ynum,                     # ynum
            xnum,                     # xnum
            ADATA.eta_flow.name,      # name
            ADATA.eta_flow.units,     # units
            ADATA.eta_flow.grid_type, # grid_type
            ADATA.eta_flow.description # description
        ),
        ScalarArray2DState(
            ynum,                      # ynum
            xnum,                      # xnum
            ADATA.eta_flow0.name,      # name
            ADATA.eta_flow0.units,     # units
            ADATA.eta_flow0.grid_type, # grid_type
            ADATA.eta_flow0.description # description
        )
    )
    update_output_format(data)
    return data
end

function update_output_format(data::Viscosity)::Nothing
    data.etas0.outform.log10 = true
    data.etas0.outform.fname = "log_etas0_Pas"
    data.etas1.outform.log10 = true
    data.etas1.outform.fname = "log_eta_Pas"
    data.etan0.outform.log10 = true
    data.etan0.outform.fname = "log_etan0_Pas"
    data.etan1.outform.log10 = true
    data.etan1.outform.fname = "log_etan1_Pas"
    data.eta_flow.outform.log10 = true
    data.eta_flow.outform.fname = "log_eta_flow_Pas"
    data.eta_flow0.outform.log10 = true
    data.eta_flow0.outform.fname = "log_eta_flow_Pas"
    return nothing
end

end # module
