module EarthLayeringGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "earth_layering"

const PDATA = get_eb_parameters()

"""
    EarthLayering <: AbstractParameterGroup

Parameter group for Earth layering geometry parameters.

# Fields
- `thick_upper_crust::`[`ParameterFloat`](@ref): $(PDATA.thick_upper_crust.description)
- `thick_lower_crust::`[`ParameterFloat`](@ref): $(PDATA.thick_lower_crust.description)
- `thick_upper_lith::`[`ParameterFloat`](@ref): $(PDATA.thick_upper_lith.description)
- `thick_middle_lith::`[`ParameterFloat`](@ref): $(PDATA.thick_middle_lith.description)
- `thick_lower_lith::`[`ParameterFloat`](@ref): $(PDATA.thick_lower_lith.description)
- `thick_lith::`[`ParameterFloat`](@ref): $(PDATA.thick_lith.description)
- `thick_crust::`[`ParameterFloat`](@ref): $(PDATA.thick_crust.description)
- `thick_mantle_lith::`[`ParameterFloat`](@ref): $(PDATA.thick_mantle_lith.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
  parameter objects

# Nested Dot Access
- `thick_upper_crust = $(ROOT_NAME).$(GRP_NAME).thick_upper_crust.value`
- `thick_lower_crust = $(ROOT_NAME).$(GRP_NAME).thick_lower_crust.value`
- `thick_upper_lith = $(ROOT_NAME).$(GRP_NAME).thick_upper_lith.value`
- `thick_middle_lith = $(ROOT_NAME).$(GRP_NAME).thick_middle_lith.value`
- `thick_lower_lith = $(ROOT_NAME).$(GRP_NAME).thick_lower_lith.value`
- `thick_lith = $(ROOT_NAME).$(GRP_NAME).thick_lith.value`
- `thick_crust = $(ROOT_NAME).$(GRP_NAME).thick_crust.value`
- `thick_mantle_lith = $(ROOT_NAME).$(GRP_NAME).thick_mantle_lith.value`

# Constructor
    EarthLayering()

Create a new EarthLayering parameter group with default values (NaN indicates unset 
parameters).

# Returns
- `EarthLayering`: New EarthLayering parameter group with initialized default values

"""
mutable struct EarthLayering <: AbstractParameterGroup
    thick_upper_crust::ParameterFloat
    thick_lower_crust::ParameterFloat
    thick_upper_lith::ParameterFloat
    thick_middle_lith::ParameterFloat
    thick_lower_lith::ParameterFloat
    thick_lith::ParameterFloat
    thick_crust::ParameterFloat
    thick_mantle_lith::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function EarthLayering()::EarthLayering
    pdata = get_eb_parameters()
    # Negative values indicate that parameter has not been set.
    data = EarthLayering(
        pdata.thick_upper_crust,
        pdata.thick_lower_crust,
        pdata.thick_upper_lith,
        pdata.thick_middle_lith,
        pdata.thick_lower_lith,
        pdata.thick_lith,
        pdata.thick_crust,
        pdata.thick_mantle_lith,
        Union{ParameterFloat}[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
