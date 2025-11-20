"""
Module for subgrid diffusion parameters.

Provides data structures for configuring numerical subgrid diffusion 
coefficients for stress and temperature interpolation.
"""
module SubgridDiffusionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.markers.parameters"
const GRP_NAME = "subgrid_diffusion"

const PDATA = get_eb_parameters()

"""
    SubgridDiffusion <: AbstractParameterGroup

Parameter group for subgrid diffusion coefficients.

# Fields
- `subgrid_diff_coef_stress::`[`ParameterFloat`](@ref): $(PDATA.subgrid_diff_coef_stress.description)
- `subgrid_diff_coef_temp::`[`ParameterFloat`](@ref): $(PDATA.subgrid_diff_coef_temp.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects.

# Nested Dot Access
- `subgrid_diff_coef_stress = $(ROOT_NAME).$(GRP_NAME).subgrid_diff_coef_stress.value`
- `subgrid_diff_coef_temp = $(ROOT_NAME).$(GRP_NAME).subgrid_diff_coef_temp.value`

# Constructor
    SubgridDiffusion()

# Returns
- A new SubgridDiffusion struct with the given marker parameters.

"""
mutable struct SubgridDiffusion <: AbstractParameterGroup
    subgrid_diff_coef_stress::ParameterFloat
    subgrid_diff_coef_temp::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function SubgridDiffusion()::SubgridDiffusion
    pdata = get_eb_parameters()
    data = SubgridDiffusion(
        pdata.subgrid_diff_coef_stress,
        pdata.subgrid_diff_coef_temp,
        ParameterFloat[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
