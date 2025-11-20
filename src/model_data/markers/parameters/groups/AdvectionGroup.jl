"""
Module for marker advection parameters.

Provides data structures for configuring marker advection schemes and 
time-stepping options.
"""
module AdvectionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.markers.parameters"
const GRP_NAME = "advection"

const PDATA = get_eb_parameters()

"""
    Advection <: AbstractParameterGroup

Parameter group for marker advection configuration.

# Fields
- `marker_cell_displ_max::`[`ParameterFloat`](@ref): $(PDATA.marker_cell_displ_max.description)
- `itype_move_markers::`[`ParameterInt`](@ref): $(PDATA.itype_move_markers.description)
- `stype_move_markers::`[`ParameterStr`](@ref): $(PDATA.stype_move_markers.description)
- `iuse_local_adaptive_time_stepping::`[`ParameterInt`](@ref): $(PDATA.iuse_local_adaptive_time_stepping.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects.

# Nested Dot Access
- `marker_cell_displ_max = $(ROOT_NAME).$(GRP_NAME).marker_cell_displ_max.value`
- `itype_move_markers = $(ROOT_NAME).$(GRP_NAME).itype_move_markers.value`
- `stype_move_markers = $(ROOT_NAME).$(GRP_NAME).stype_move_markers.value`
- `iuse_local_adaptive_time_stepping = $(ROOT_NAME).$(GRP_NAME).iuse_local_adaptive_time_stepping.value`

# Constructor
    Advection()

# Returns
- A new Advection struct with the given marker parameters.

"""
mutable struct Advection <: AbstractParameterGroup
    marker_cell_displ_max::ParameterFloat
    itype_move_markers::ParameterInt
    stype_move_markers::ParameterStr
    iuse_local_adaptive_time_stepping::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Advection()::Advection
    pdata = get_eb_parameters()
    data = Advection(
        pdata.marker_cell_displ_max,
        pdata.itype_move_markers,
        pdata.stype_move_markers,
        pdata.iuse_local_adaptive_time_stepping,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
