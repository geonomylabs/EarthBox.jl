"""
Module for marker recycling parameters.

Provides data structures for configuring marker recycling at domain boundaries.
"""
module RecyclingGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.markers.parameters"
const GRP_NAME = "recycling"

const PDATA = get_eb_parameters()

"""
    Recycling <: AbstractParameterGroup

Parameter group for marker recycling configuration.

# Fields
- `imantle::`[`ParameterInt`](@ref): $(PDATA.imantle.description)
- `obj_list::Vector{ParameterInt}`: List of numerical parameter objects

# Nested Dot Access
- `imantle = $(ROOT_NAME).$(GRP_NAME).imantle.value`

# Constructor
    Recycling()

# Returns
- A new Recycling struct with the given marker parameters.

"""
mutable struct Recycling <: AbstractParameterGroup
    imantle::ParameterInt
    obj_list::Vector{ParameterInt}
end

function Recycling()::Recycling
    pdata = get_eb_parameters()
    data = Recycling(
        pdata.imantle,
        ParameterInt[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
