module VelocityCalcOptionsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.stokes_continuity.parameters"
const GRP_NAME = "velocity_calc_option"

const PDATA = get_eb_parameters()

"""
    VelocityCalcOptions <: AbstractParameterGroup

Parameter group for velocity calculation options.

# Fields
- `itype_velocity::`[`ParameterInt`](@ref): $(PDATA.itype_velocity.description)
- `stype_velocity::`[`ParameterStr`](@ref): $(PDATA.stype_velocity.description)
- `obj_list::Vector{ParameterInt}`: List of parameter objects

# Nested Dot Access
- `itype_velocity = $(ROOT_NAME).$(GRP_NAME).itype_velocity.value`
- `stype_velocity = $(ROOT_NAME).$(GRP_NAME).stype_velocity.value`

# Constructor
    VelocityCalcOptions()

# Returns
- `VelocityCalcOptions`: New VelocityCalcOptions parameter group with initialized values

"""
mutable struct VelocityCalcOptions <: AbstractParameterGroup
    itype_velocity::ParameterInt
    stype_velocity::ParameterStr
    obj_list::Vector{ParameterInt}
end

function VelocityCalcOptions()::VelocityCalcOptions
    pdata = get_eb_parameters()
    data = VelocityCalcOptions(
        pdata.itype_velocity,
        pdata.stype_velocity,
        ParameterInt[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
