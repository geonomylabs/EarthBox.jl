"""
Module for multiple time step increase parameters.

Provides data structures for configuring multiple increases of the 
viscoelastic time step during a simulation at specified time step numbers.
"""
module MultipleIncreaseGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.timestep.parameters"
const GRP_NAME = "multiple_increase"

const PDATA = get_eb_parameters()

"""
    MultipleIncrease <: AbstractParameterGroup

Parameter group for multiple time step increases.

# Fields
- `iuse_multiple_timestep_increase::`[`ParameterInt`](@ref): $(PDATA.iuse_multiple_timestep_increase.description)
- `ntime_increase_1::`[`ParameterInt`](@ref): $(PDATA.ntime_increase_1.description)
- `ntime_increase_2::`[`ParameterInt`](@ref): $(PDATA.ntime_increase_2.description)
- `ntime_increase_3::`[`ParameterInt`](@ref): $(PDATA.ntime_increase_3.description)
- `ntime_increase_4::`[`ParameterInt`](@ref): $(PDATA.ntime_increase_4.description)
- `time_increase_factor::`[`ParameterFloat`](@ref): $(PDATA.time_increase_factor.description)
- `cell_displ_factor::`[`ParameterFloat`](@ref): $(PDATA.cell_displ_factor.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of numerical parameter objects

# Nested Dot Access
- `iuse_multiple_timestep_increase = $(ROOT_NAME).$(GRP_NAME).iuse_multiple_timestep_increase.value`
- `ntime_increase_1 = $(ROOT_NAME).$(GRP_NAME).ntime_increase_1.value`
- `ntime_increase_2 = $(ROOT_NAME).$(GRP_NAME).ntime_increase_2.value`
- `ntime_increase_3 = $(ROOT_NAME).$(GRP_NAME).ntime_increase_3.value`
- `ntime_increase_4 = $(ROOT_NAME).$(GRP_NAME).ntime_increase_4.value`
- `time_increase_factor = $(ROOT_NAME).$(GRP_NAME).time_increase_factor.value`
- `cell_displ_factor = $(ROOT_NAME).$(GRP_NAME).cell_displ_factor.value`

# Constructor
    MultipleIncrease()

Initializes multiple time step increase parameters. Time increase factor and 
cell displacement factor are both set to 2.0.
"""
mutable struct MultipleIncrease <: AbstractParameterGroup
    iuse_multiple_timestep_increase::ParameterInt
    ntime_increase_1::ParameterInt
    ntime_increase_2::ParameterInt
    ntime_increase_3::ParameterInt
    ntime_increase_4::ParameterInt
    time_increase_factor::ParameterFloat
    cell_displ_factor::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function MultipleIncrease()::MultipleIncrease
    pdata = get_eb_parameters()
    data = MultipleIncrease(
        pdata.iuse_multiple_timestep_increase,
        pdata.ntime_increase_1,
        pdata.ntime_increase_2,
        pdata.ntime_increase_3,
        pdata.ntime_increase_4,
        pdata.time_increase_factor,
        pdata.cell_displ_factor,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
