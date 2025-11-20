module MobileWallGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "mobile_wall"

const PDATA = get_eb_parameters()

"""
    MobileWall <: AbstractParameterGroup

Parameter group for mobile wall geometry parameters.

# Fields
- `x_left_mobile_wall::`[`ParameterFloat`](@ref): $(PDATA.x_left_mobile_wall.description)
- `x_right_mobile_wall::`[`ParameterFloat`](@ref): $(PDATA.x_right_mobile_wall.description)
- `y_top_mobile_wall::`[`ParameterFloat`](@ref): $(PDATA.y_top_mobile_wall.description)
- `y_bottom_mobile_wall::`[`ParameterFloat`](@ref): $(PDATA.y_bottom_mobile_wall.description)
- `plate_extension_width::`[`ParameterFloat`](@ref): $(PDATA.plate_extension_width.description)
- `plate_extension_thickness::`[`ParameterFloat`](@ref): $(PDATA.plate_extension_thickness.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `x_left_mobile_wall = $(ROOT_NAME).$(GRP_NAME).x_left_mobile_wall.value`
- `x_right_mobile_wall = $(ROOT_NAME).$(GRP_NAME).x_right_mobile_wall.value`
- `y_top_mobile_wall = $(ROOT_NAME).$(GRP_NAME).y_top_mobile_wall.value`
- `y_bottom_mobile_wall = $(ROOT_NAME).$(GRP_NAME).y_bottom_mobile_wall.value`
- `plate_extension_width = $(ROOT_NAME).$(GRP_NAME).plate_extension_width.value`
- `plate_extension_thickness = $(ROOT_NAME).$(GRP_NAME).plate_extension_thickness.value`

# Constructor
    MobileWall()

Create a new MobileWall parameter group with default values.

# Returns
- `MobileWall`: New MobileWall parameter group with initialized default values

"""
mutable struct MobileWall <: AbstractParameterGroup
    x_left_mobile_wall::ParameterFloat
    x_right_mobile_wall::ParameterFloat
    y_top_mobile_wall::ParameterFloat
    y_bottom_mobile_wall::ParameterFloat
    plate_extension_width::ParameterFloat
    plate_extension_thickness::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function MobileWall()::MobileWall
    pdata = get_eb_parameters()
    data = MobileWall(
        pdata.x_left_mobile_wall,
        pdata.x_right_mobile_wall,
        pdata.y_top_mobile_wall,
        pdata.y_bottom_mobile_wall,
        pdata.plate_extension_width,
        pdata.plate_extension_thickness,
        ParameterFloat[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
