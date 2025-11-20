module StickyAirGeometryGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.geometry.parameters"
const GRP_NAME = "sticky_air_geometry"

const PDATA = get_eb_parameters()

"""
    StickyAirGeometry <: AbstractParameterGroup

Parameter group for sticky air geometry parameters.

# Fields
- `thick_air::`[`ParameterFloat`](@ref): $(PDATA.thick_air.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `thick_air = $(ROOT_NAME).$(GRP_NAME).thick_air.value`

# Constructor
    StickyAirGeometry()

Create a new StickyAirGeometry parameter group with default values.

# Returns
- `StickyAirGeometry`: New StickyAirGeometry parameter group with initialized default values

"""
mutable struct StickyAirGeometry <: AbstractParameterGroup
    thick_air::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function StickyAirGeometry()::StickyAirGeometry
    pdata = get_eb_parameters()
    data = StickyAirGeometry(
        pdata.thick_air,
        ParameterFloat[]
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
