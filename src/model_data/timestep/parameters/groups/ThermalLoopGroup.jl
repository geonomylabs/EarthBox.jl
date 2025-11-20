"""
Module for thermal loop parameters.

Provides data structures for configuring the adaptive heat solver loop 
including thermal time steps and time tracking.
"""
module ThermalLoopGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.timestep.parameters"
const GRP_NAME = "thermal_loop"

const PDATA = get_eb_parameters()

"""
    ThermalLoop <: AbstractParameterGroup

Parameter group for thermal loop configuration.

# Fields
- `timestep_heat::`[`ParameterFloat`](@ref): $(PDATA.timestep_heat.description)
- `timestep_sum::`[`ParameterFloat`](@ref): $(PDATA.timestep_sum.description)
- `obj_list::Vector{ParameterFloat}`: List of numerical parameter objects

# Nested Dot Access
- `timestep_heat = $(ROOT_NAME).$(GRP_NAME).timestep_heat.value`
- `timestep_sum = $(ROOT_NAME).$(GRP_NAME).timestep_sum.value`

# Constructor
    ThermalLoop()

Initializes thermal loop parameters with zero values for thermal time step 
and total time.
"""
mutable struct ThermalLoop <: AbstractParameterGroup
    timestep_heat::ParameterFloat
    timestep_sum::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function ThermalLoop()::ThermalLoop
    pdata = get_eb_parameters()
    data = ThermalLoop(
        pdata.timestep_heat,
        pdata.timestep_sum,
        ParameterFloat[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
