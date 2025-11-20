module OutputStokesGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.stokes_continuity.parameters"
const GRP_NAME = "output"

const PDATA = get_eb_parameters()

"""
    OutputStokes <: AbstractParameterGroup

Parameter group for Stokes solver output parameters.

# Fields
- `outtest::`[`ParameterInt`](@ref): $(PDATA.outtest.description)
- `outtest2::`[`ParameterInt`](@ref): $(PDATA.outtest2.description)
- `outtest3::`[`ParameterInt`](@ref): $(PDATA.outtest3.description)
- `obj_list::Vector{ParameterInt}`: List of parameter objects

# Nested Dot Access
- `outtest = $(ROOT_NAME).$(GRP_NAME).outtest.value`
- `outtest2 = $(ROOT_NAME).$(GRP_NAME).outtest2.value`
- `outtest3 = $(ROOT_NAME).$(GRP_NAME).outtest3.value`

# Constructor
    OutputStokes()

# Returns
- `OutputStokes`: New OutputStokes parameter group with initialized values

"""
mutable struct OutputStokes <: AbstractParameterGroup
    outtest::ParameterInt
    outtest2::ParameterInt
    outtest3::ParameterInt
    obj_list::Vector{ParameterInt}
end

function OutputStokes()::OutputStokes
    pdata = get_eb_parameters()
    data = OutputStokes(
        pdata.outtest,
        pdata.outtest2,
        pdata.outtest3,
        ParameterInt[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
