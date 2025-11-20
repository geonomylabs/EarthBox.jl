module RheologyGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.melting.parameters"
const GRP_NAME = "rheology"

const PDATA = get_eb_parameters()

"""
    Rheology <: AbstractParameterGroup

Parameter group for rheological properties of molten rock.

# Fields
- `viscosity_melt::`[`ParameterFloat`](@ref): $(PDATA.viscosity_melt.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `viscosity_melt = $(ROOT_NAME).$(GRP_NAME).viscosity_melt.value`

# Constructor
    Rheology()

# Returns
- `Rheology`: New Rheology parameter group with initialized values

"""
mutable struct Rheology <: AbstractParameterGroup
    viscosity_melt::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function Rheology()::Rheology
    pdata = get_eb_parameters()
    data = Rheology(
        pdata.viscosity_melt,
        ParameterFloat[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
