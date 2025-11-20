"""
Module for pressure boundary condition parameters.

Provides data structures for configuring pressure boundary conditions.
"""
module PressureGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: AbstractParameterGroup
import EarthBox.Parameters: ParameterFloat

const ROOT_NAME = "model.bcs.parameters"
const GRP_NAME = "pressure"

const PDATA = get_eb_parameters()

"""
    Pressure <: AbstractParameterGroup

Parameter group for pressure boundary conditions.

# Fields
- `pressure_bc::`[`ParameterFloat`](@ref): $(PDATA.pressure_bc.description)

# Nested Dot Access
- `pressure_bc = $(ROOT_NAME).$(GRP_NAME).pressure_bc.value`

# Constructor
    Pressure()

Initializes pressure boundary condition parameter with default value of 
10000.0 Pa.
"""
mutable struct Pressure <: AbstractParameterGroup
    pressure_bc::ParameterFloat
end

function Pressure()::Pressure
    pdata = get_eb_parameters()
    data = Pressure(
        pdata.pressure_bc
    )
    return data
end

end # module
