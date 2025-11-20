"""
Module for main time loop parameters.

Provides data structures for configuring the main time-stepping loop including
model duration, time step sizes, and viscoelastic time stepping.
"""
module MainTimeLoopGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.timestep.parameters"
const GRP_NAME = "main_time_loop"

const PDATA = get_eb_parameters()

"""
    MainTimeLoop <: AbstractParameterGroup

Parameter group for main time loop configuration.

# Fields
- `model_duration_myr::`[`ParameterFloat`](@ref): $(PDATA.model_duration_myr.description)
- `ntimestep_max::`[`ParameterInt`](@ref): $(PDATA.ntimestep_max.description)
- `timesum::`[`ParameterFloat`](@ref): $(PDATA.timesum.description)
- `ntimestep::`[`ParameterInt`](@ref): $(PDATA.ntimestep.description)
- `timestep_viscoelastic::`[`ParameterFloat`](@ref): $(PDATA.timestep_viscoelastic.description)
- `timestep_viscoelastic_step1::`[`ParameterFloat`](@ref): $(PDATA.timestep_viscoelastic_step1.description)
- `timestep_viscoelastic_step2::`[`ParameterFloat`](@ref): $(PDATA.timestep_viscoelastic_step2.description)
- `timestep::`[`ParameterFloat`](@ref): $(PDATA.timestep.description)
- `iupdate_timestep::`[`ParameterInt`](@ref): $(PDATA.iupdate_timestep.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `model_duration_myr = $(ROOT_NAME).$(GRP_NAME).model_duration_myr.value`
- `ntimestep_max = $(ROOT_NAME).$(GRP_NAME).ntimestep_max.value`
- `timesum = $(ROOT_NAME).$(GRP_NAME).timesum.value`
- `ntimestep = $(ROOT_NAME).$(GRP_NAME).ntimestep.value`
- `timestep_viscoelastic = $(ROOT_NAME).$(GRP_NAME).timestep_viscoelastic.value`
- `timestep_viscoelastic_step1 = $(ROOT_NAME).$(GRP_NAME).timestep_viscoelastic_step1.value`
- `timestep_viscoelastic_step2 = $(ROOT_NAME).$(GRP_NAME).timestep_viscoelastic_step2.value`
- `timestep = $(ROOT_NAME).$(GRP_NAME).timestep.value`
- `iupdate_timestep = $(ROOT_NAME).$(GRP_NAME).iupdate_timestep.value`

# Constructor
    MainTimeLoop()

Initializes main time loop parameters with default values. Main time step is 
set to 10000.0 s.
"""
mutable struct MainTimeLoop <: AbstractParameterGroup
    model_duration_myr::ParameterFloat
    ntimestep_max::ParameterInt
    timesum::ParameterFloat
    ntimestep::ParameterInt
    timestep_viscoelastic::ParameterFloat
    timestep_viscoelastic_step1::ParameterFloat
    timestep_viscoelastic_step2::ParameterFloat
    timestep::ParameterFloat
    iupdate_timestep::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function MainTimeLoop()::MainTimeLoop
    pdata = get_eb_parameters()
    data = MainTimeLoop(
        pdata.model_duration_myr,
        pdata.ntimestep_max,
        pdata.timesum,
        pdata.ntimestep,
        pdata.timestep_viscoelastic,
        pdata.timestep_viscoelastic_step1,
        pdata.timestep_viscoelastic_step2,
        pdata.timestep,
        pdata.iupdate_timestep,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
