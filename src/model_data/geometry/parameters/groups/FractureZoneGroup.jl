module FractureZoneGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "fracture_zone"

const PDATA = get_eb_parameters()

"""
    FractureZone <: AbstractParameterGroup

Parameter group for fracture zone geometry parameters.

# Fields
- `sediment_thickness::`[`ParameterFloat`](@ref): $(PDATA.sediment_thickness.description)
- `basaltic_oceanic_crust_thickness::`[`ParameterFloat`](@ref): $(PDATA.basaltic_oceanic_crust_thickness.description)
- `gabbroic_oceanic_crust_thickness::`[`ParameterFloat`](@ref): $(PDATA.gabbroic_oceanic_crust_thickness.description)
- `thickness_of_younger_lithosphere::`[`ParameterFloat`](@ref): $(PDATA.thickness_of_younger_lithosphere.description)
- `thickness_of_older_lithosphere::`[`ParameterFloat`](@ref): $(PDATA.thickness_of_older_lithosphere.description)
- `thickness_of_weak_lithosphere::`[`ParameterFloat`](@ref): $(PDATA.thickness_of_weak_lithosphere.description)
- `x_fracture_zone_start::`[`ParameterFloat`](@ref): $(PDATA.x_fracture_zone_start.description)
- `x_fracture_zone_end::`[`ParameterFloat`](@ref): $(PDATA.x_fracture_zone_end.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
  parameter objects

# Nested Dot Access
- `sediment_thickness = $(ROOT_NAME).$(GRP_NAME).sediment_thickness.value`
- `basaltic_oceanic_crust_thickness = 
    $(ROOT_NAME).$(GRP_NAME).basaltic_oceanic_crust_thickness.value`
- `gabbroic_oceanic_crust_thickness = 
    $(ROOT_NAME).$(GRP_NAME).gabbroic_oceanic_crust_thickness.value`
- `thickness_of_younger_lithosphere = 
  $(ROOT_NAME).$(GRP_NAME).thickness_of_younger_lithosphere.value`
- `thickness_of_older_lithosphere = 
   $(ROOT_NAME).$(GRP_NAME).thickness_of_older_lithosphere.value`
- `thickness_of_weak_lithosphere = 
   $(ROOT_NAME).$(GRP_NAME).thickness_of_weak_lithosphere.value`
- `x_fracture_zone_start = 
   $(ROOT_NAME).$(GRP_NAME).x_fracture_zone_start.value`
- `x_fracture_zone_end = $(ROOT_NAME).$(GRP_NAME).x_fracture_zone_end.value`

# Constructor
    FractureZone()

Create a new FractureZone parameter group with default values.

# Returns
- `FractureZone`: New FractureZone parameter group with initialized default values

"""
mutable struct FractureZone <: AbstractParameterGroup
    sediment_thickness::ParameterFloat
    basaltic_oceanic_crust_thickness::ParameterFloat
    gabbroic_oceanic_crust_thickness::ParameterFloat
    thickness_of_younger_lithosphere::ParameterFloat
    thickness_of_older_lithosphere::ParameterFloat
    thickness_of_weak_lithosphere::ParameterFloat
    x_fracture_zone_start::ParameterFloat
    x_fracture_zone_end::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function FractureZone()::FractureZone
    pdata = get_eb_parameters()
    data = FractureZone(
        pdata.sediment_thickness,
        pdata.basaltic_oceanic_crust_thickness,
        pdata.gabbroic_oceanic_crust_thickness,
        pdata.thickness_of_younger_lithosphere,
        pdata.thickness_of_older_lithosphere,
        pdata.thickness_of_weak_lithosphere,
        pdata.x_fracture_zone_start,
        pdata.x_fracture_zone_end,
        Union{ParameterFloat}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 