module WeakFaultGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "weak_fault"

const PDATA = get_eb_parameters()

"""
    WeakFault <: AbstractParameterGroup

Parameter group for weak fault geometry parameters.

# Fields
- `fault_dip_degrees::`[`ParameterFloat`](@ref): $(PDATA.fault_dip_degrees.description)
- `fault_thickness::`[`ParameterFloat`](@ref): $(PDATA.fault_thickness.description)
- `x_initial_fault::`[`ParameterFloat`](@ref): $(PDATA.x_initial_fault.description)
- `fault_height::`[`ParameterFloat`](@ref): $(PDATA.fault_height.description)
- `iuse_weak_fault::`[`ParameterInt`](@ref): $(PDATA.iuse_weak_fault.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
    parameter objects

# Nested Dot Access
- `fault_dip_degrees = $(ROOT_NAME).$(GRP_NAME).fault_dip_degrees.value`
- `fault_thickness = $(ROOT_NAME).$(GRP_NAME).fault_thickness.value`
- `x_initial_fault = $(ROOT_NAME).$(GRP_NAME).x_initial_fault.value`
- `fault_height = $(ROOT_NAME).$(GRP_NAME).fault_height.value`
- `iuse_weak_fault = $(ROOT_NAME).$(GRP_NAME).iuse_weak_fault.value`

# Constructor
    WeakFault()

Create a new WeakFault parameter group with default values.

# Returns
- `WeakFault`: New WeakFault parameter group with initialized default values

"""
mutable struct WeakFault <: AbstractParameterGroup
    fault_dip_degrees::ParameterFloat
    fault_thickness::ParameterFloat
    x_initial_fault::ParameterFloat
    fault_height::ParameterFloat
    iuse_weak_fault::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function WeakFault()::WeakFault
    pdata = get_eb_parameters()
    data = WeakFault(
        pdata.fault_dip_degrees,
        pdata.fault_thickness,
        pdata.x_initial_fault,
        pdata.fault_height,
        pdata.iuse_weak_fault,
        Union{ParameterFloat, ParameterInt}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 