"""
Module for output step parameters.

Provides data structures for configuring output frequency and tracking output
generation during model runs.
"""
module OutputStepsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.timestep.parameters"
const GRP_NAME = "output_steps"

const PDATA = get_eb_parameters()

"""
    OutputSteps <: AbstractParameterGroup

Parameter group for output step configuration.

# Fields
- `timestep_out::`[`ParameterFloat`](@ref): $(PDATA.timestep_out.description)
- `nskip::`[`ParameterInt`](@ref): $(PDATA.nskip.description)
- `icount_output::`[`ParameterInt`](@ref): $(PDATA.icount_output.description)
- `noutput::`[`ParameterInt`](@ref): $(PDATA.noutput.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of numerical parameter objects

# Nested Dot Access
- `timestep_out = $(ROOT_NAME).$(GRP_NAME).timestep_out.value`
- `nskip = $(ROOT_NAME).$(GRP_NAME).nskip.value`
- `icount_output = $(ROOT_NAME).$(GRP_NAME).icount_output.value`
- `noutput = $(ROOT_NAME).$(GRP_NAME).noutput.value`

# Constructor
    OutputSteps()

Initializes output step parameters with default values. Output time step is 
set to 1000.0 s.
"""
mutable struct OutputSteps <: AbstractParameterGroup
    timestep_out::ParameterFloat
    nskip::ParameterInt
    icount_output::ParameterInt
    noutput::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function OutputSteps()::OutputSteps
    pdata = get_eb_parameters()
    data = OutputSteps(
        pdata.timestep_out,
        pdata.nskip,
        pdata.icount_output,
        pdata.noutput,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
