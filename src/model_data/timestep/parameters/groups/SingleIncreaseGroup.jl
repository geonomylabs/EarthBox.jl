"""
Module for single time step increase parameters.

Provides data structures for configuring a single increase of the time step
during a simulation at a specified time step number.
"""
module SingleIncreaseGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.timestep.parameters"
const GRP_NAME = "single_increase"

const PDATA = get_eb_parameters()

"""
    SingleIncrease <: AbstractParameterGroup

Parameter group for single time step increase.

# Fields
- `iuse_single_timestep_increase::`[`ParameterInt`](@ref): $(PDATA.iuse_single_timestep_increase.description)
- `ntime_increase::`[`ParameterInt`](@ref): $(PDATA.ntime_increase.description)
- `obj_list::Vector{ParameterInt}`: List of parameter objects

# Nested Dot Access
- `iuse_single_timestep_increase = $(ROOT_NAME).$(GRP_NAME).iuse_single_timestep_increase.value`
- `ntime_increase = $(ROOT_NAME).$(GRP_NAME).ntime_increase.value`

# Constructor
    SingleIncrease()

Initializes single time step increase parameters. The increase is set to 
occur at time step 100.
"""
mutable struct SingleIncrease <: AbstractParameterGroup
    iuse_single_timestep_increase::ParameterInt
    ntime_increase::ParameterInt
    obj_list::Vector{ParameterInt}
end

function SingleIncrease()::SingleIncrease
    pdata = get_eb_parameters()
    data = SingleIncrease(
        pdata.iuse_single_timestep_increase,
        pdata.ntime_increase,
        ParameterInt[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
