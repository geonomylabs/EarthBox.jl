module BoundaryFrictionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "boundary_friction"

const PDATA = get_eb_parameters()

"""
    BoundaryFriction <: AbstractParameterGroup

Parameter group for boundary friction zone properties.

# Fields
- `boundary_friction_width::`[`ParameterFloat`](@ref): $(PDATA.boundary_friction_width.description)
- `boundary_friction_angle::`[`ParameterFloat`](@ref): $(PDATA.boundary_friction_angle.description)
- `boundary_cohesion::`[`ParameterFloat`](@ref): $(PDATA.boundary_cohesion.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `boundary_friction_width = $(ROOT_NAME).$(GRP_NAME).boundary_friction_width.value`
- `boundary_friction_angle = $(ROOT_NAME).$(GRP_NAME).boundary_friction_angle.value`
- `boundary_cohesion = $(ROOT_NAME).$(GRP_NAME).boundary_cohesion.value`

# Constructor
    BoundaryFriction()

Create a new BoundaryFriction parameter group with default values.

# Returns
- `BoundaryFriction`: New BoundaryFriction parameter group with initialized values

"""
mutable struct BoundaryFriction <: AbstractParameterGroup
    boundary_friction_width::ParameterFloat
    boundary_friction_angle::ParameterFloat
    boundary_cohesion::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function BoundaryFriction()::BoundaryFriction
    pdata = get_eb_parameters()
    data = BoundaryFriction(
        pdata.boundary_friction_width,
        pdata.boundary_friction_angle,
        pdata.boundary_cohesion,
        ParameterFloat[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
