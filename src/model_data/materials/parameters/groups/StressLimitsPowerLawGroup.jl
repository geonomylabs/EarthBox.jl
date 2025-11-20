module StressLimitsPowerLawGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "stress_limits_powerlaw"

const PDATA = get_eb_parameters()

"""
    StressLimitsPowerLaw <: AbstractParameterGroup

Parameter group for power-law stress limits.

# Fields
- `powerlaw_stress_min::`[`ParameterFloat`](@ref): $(PDATA.powerlaw_stress_min.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `powerlaw_stress_min = $(ROOT_NAME).$(GRP_NAME).powerlaw_stress_min.value`

# Constructor
    StressLimitsPowerLaw()

# Returns
- `StressLimitsPowerLaw`: New StressLimitsPowerLaw parameter group with initialized values

"""
mutable struct StressLimitsPowerLaw <: AbstractParameterGroup
    powerlaw_stress_min::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function StressLimitsPowerLaw()::StressLimitsPowerLaw
    pdata = get_eb_parameters()
    data = StressLimitsPowerLaw(
        pdata.powerlaw_stress_min,
        ParameterFloat[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
