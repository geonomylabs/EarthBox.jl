"""
Module for boundary condition option parameters.

Provides data structures for configuring boundary condition options for Stokes
and heat equations.
"""
module BcOptionsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: AbstractParameterGroup
import EarthBox.Parameters: ParameterInt, ParameterStr, ParameterFloat
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list

const ROOT_NAME = "model.bcs.parameters"
const GRP_NAME = "bc_options"

const PDATA = get_eb_parameters()

"""
    BcOptions <: AbstractParameterGroup

Parameter group for boundary condition options.

# Fields
- `itype_bc::`[`ParameterInt`](@ref): $(PDATA.itype_bc.description)
- `stype_bc::`[`ParameterStr`](@ref): $(PDATA.stype_bc.description)
- `pressure_bc_mode::`[`ParameterFloat`](@ref): $(PDATA.pressure_bc_mode.description)

# Nested Dot Access
- `itype_bc = $(ROOT_NAME).$(GRP_NAME).itype_bc.value`
- `stype_bc = $(ROOT_NAME).$(GRP_NAME).stype_bc.value`
- `pressure_bc_mode = $(ROOT_NAME).$(GRP_NAME).pressure_bc_mode.value`

# Constructor
    BcOptions()

Initializes boundary condition option parameters with default values.
"""
mutable struct BcOptions <: AbstractParameterGroup
    itype_bc::ParameterInt
    stype_bc::ParameterStr
    pressure_bc_mode::ParameterFloat
end

function BcOptions()::BcOptions
    pdata = get_eb_parameters()
    data = BcOptions(
        pdata.itype_bc,
        pdata.stype_bc,
        pdata.pressure_bc_mode
    )
    return data
end

end # module
