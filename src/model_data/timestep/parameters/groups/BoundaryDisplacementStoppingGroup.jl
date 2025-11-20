"""
Module for boundary displacement stopping criteria parameters.

Provides data structures for configuring model stopping conditions based on
boundary displacement or extensional strain limits.
"""
module BoundaryDisplacementStoppingGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.timestep.parameters"
const GRP_NAME = "boundary_displacement_stopping"

const PDATA = get_eb_parameters()

"""
    BoundaryDisplacementStopping <: AbstractParameterGroup

Parameter group for stopping criteria for boundary displacement or extensional strain.

# Fields
- `iuse_boundary_displacement::`[`ParameterInt`](@ref): $(PDATA.iuse_boundary_displacement.description)
- `displ_limit::`[`ParameterFloat`](@ref): $(PDATA.displ_limit.description)
- `iuse_extensional_strain::`[`ParameterInt`](@ref): $(PDATA.iuse_extensional_strain.description)
- `strain_limit::`[`ParameterFloat`](@ref): $(PDATA.strain_limit.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of numerical 
    parameter objects

# Nested Dot Access
- `iuse_boundary_displacement = $(ROOT_NAME).$(GRP_NAME).iuse_boundary_displacement.value`
- `displ_limit = $(ROOT_NAME).$(GRP_NAME).displ_limit.value`
- `iuse_extensional_strain = $(ROOT_NAME).$(GRP_NAME).iuse_extensional_strain.value`
- `strain_limit = $(ROOT_NAME).$(GRP_NAME).strain_limit.value`

# Constructor
    BoundaryDisplacementStopping()

"""
mutable struct BoundaryDisplacementStopping <: AbstractParameterGroup
    iuse_boundary_displacement::ParameterInt
    displ_limit::ParameterFloat
    iuse_extensional_strain::ParameterInt
    strain_limit::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function BoundaryDisplacementStopping()::BoundaryDisplacementStopping
    pdata = get_eb_parameters()
    data = BoundaryDisplacementStopping(
        pdata.iuse_boundary_displacement,
        pdata.displ_limit,
        pdata.iuse_extensional_strain,
        pdata.strain_limit,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
