module SandboxGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "sandbox"

const PDATA = get_eb_parameters()

"""
    Sandbox <: AbstractParameterGroup

Parameter group for sandbox model geometry parameters.

# Fields
- `nsand_layers::`[`ParameterInt`](@ref): $(PDATA.nsand_layers.description)
- `y_sand_air_interface::`[`ParameterFloat`](@ref): $(PDATA.y_sand_air_interface.description)
- `y_top_microbeads::`[`ParameterFloat`](@ref): $(PDATA.y_top_microbeads.description)
- `y_bottom_microbeads::`[`ParameterFloat`](@ref): $(PDATA.y_bottom_microbeads.description)
- `x_left_ramp::`[`ParameterFloat`](@ref): $(PDATA.x_left_ramp.description)
- `x_right_ramp::`[`ParameterFloat`](@ref): $(PDATA.x_right_ramp.description)
- `ramp_dip_deg::`[`ParameterFloat`](@ref): $(PDATA.ramp_dip_deg.description)
- `pdms_layer_width::`[`ParameterFloat`](@ref): $(PDATA.pdms_layer_width.description)
- `pdms_layer_thickness::`[`ParameterFloat`](@ref): $(PDATA.pdms_layer_thickness.description)
- `obj_list::Vector{Union{`[`ParameterFloat`](@ref), [`ParameterInt`](@ref)}}`: List of 
    parameter objects

# Nested Dot Access
- `nsand_layers = $(ROOT_NAME).$(GRP_NAME).nsand_layers.value`
- `y_sand_air_interface = $(ROOT_NAME).$(GRP_NAME).y_sand_air_interface.value`
- `y_top_microbeads = $(ROOT_NAME).$(GRP_NAME).y_top_microbeads.value`
- `y_bottom_microbeads = $(ROOT_NAME).$(GRP_NAME).y_bottom_microbeads.value`
- `x_left_ramp = $(ROOT_NAME).$(GRP_NAME).x_left_ramp.value`
- `x_right_ramp = $(ROOT_NAME).$(GRP_NAME).x_right_ramp.value`
- `ramp_dip_deg = $(ROOT_NAME).$(GRP_NAME).ramp_dip_deg.value`
- `pdms_layer_width = $(ROOT_NAME).$(GRP_NAME).pdms_layer_width.value`
- `pdms_layer_thickness = $(ROOT_NAME).$(GRP_NAME).pdms_layer_thickness.value`

# Constructor
    Sandbox()

Create a new Sandbox parameter group with default values.

# Returns
- `Sandbox`: New Sandbox parameter group with initialized default values

"""
mutable struct Sandbox <: AbstractParameterGroup
    nsand_layers::ParameterInt
    y_sand_air_interface::ParameterFloat
    y_top_microbeads::ParameterFloat
    y_bottom_microbeads::ParameterFloat
    x_left_ramp::ParameterFloat
    x_right_ramp::ParameterFloat
    ramp_dip_deg::ParameterFloat
    pdms_layer_width::ParameterFloat
    pdms_layer_thickness::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Sandbox()::Sandbox
    pdata = get_eb_parameters()
    data = Sandbox(
        pdata.nsand_layers,
        pdata.y_sand_air_interface,
        pdata.y_top_microbeads,
        pdata.y_bottom_microbeads,
        pdata.x_left_ramp,
        pdata.x_right_ramp,
        pdata.ramp_dip_deg,
        pdata.pdms_layer_width,
        pdata.pdms_layer_thickness,
        Union{ParameterFloat, ParameterInt}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
