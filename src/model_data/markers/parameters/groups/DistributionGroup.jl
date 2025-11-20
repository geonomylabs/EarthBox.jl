"""
Module for marker distribution parameters.

Provides data structures for configuring the initial spatial distribution of
markers including spacing, density, and randomization options.
"""
module DistributionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import ....Grids2dContainer: Grids
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.markers.parameters"
const GRP_NAME = "distribution"

const PDATA = get_eb_parameters()

"""
    Distribution <: AbstractParameterGroup

Parameter group for marker distribution configuration.

# Fields
- `iuse_random::`[`ParameterInt`](@ref): $(PDATA.iuse_random.description)
- `dx_marker::`[`ParameterFloat`](@ref): $(PDATA.dx_marker.description)
- `dy_marker::`[`ParameterFloat`](@ref): $(PDATA.dy_marker.description)
- `nmarkers_cell_x::`[`ParameterFloat`](@ref): $(PDATA.nmarkers_cell_x.description)
- `nmarkers_cell_y::`[`ParameterFloat`](@ref): $(PDATA.nmarkers_cell_y.description)
- `mxnum::`[`ParameterInt`](@ref): $(PDATA.mxnum.description)
- `mynum::`[`ParameterInt`](@ref): $(PDATA.mynum.description)
- `marknum::`[`ParameterInt`](@ref): $(PDATA.marknum.description)
- `mxstep::`[`ParameterFloat`](@ref): $(PDATA.mxstep.description)
- `mystep::`[`ParameterFloat`](@ref): $(PDATA.mystep.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects.

# Nested Dot Access
- `iuse_random = $(ROOT_NAME).$(GRP_NAME).iuse_random.value`
- `dx_marker = $(ROOT_NAME).$(GRP_NAME).dx_marker.value`
- `dy_marker = $(ROOT_NAME).$(GRP_NAME).dy_marker.value`
- `nmarkers_cell_x = $(ROOT_NAME).$(GRP_NAME).nmarkers_cell_x.value`
- `nmarkers_cell_y = $(ROOT_NAME).$(GRP_NAME).nmarkers_cell_y.value`
- `mxnum = $(ROOT_NAME).$(GRP_NAME).mxnum.value`
- `mynum = $(ROOT_NAME).$(GRP_NAME).mynum.value`
- `marknum = $(ROOT_NAME).$(GRP_NAME).marknum.value`
- `mxstep = $(ROOT_NAME).$(GRP_NAME).mxstep.value`
- `mystep = $(ROOT_NAME).$(GRP_NAME).mystep.value`

# Constructor
    Distribution(marker_parameters::NamedTuple)

Initializes marker distribution parameters from the provided named tuple.

## Arguments
- `marker_parameters::NamedTuple`: Named tuple containing marker distribution 
    parameters.

The named tuple `marker_parameters` must contain the following parameters:
- `dx_marker`: Initial average marker spacing in horizontal direction in meters.
- `dy_marker`: Initial average marker spacing in vertical direction in meters.
- `nmarkers_cell_x`: Number of markers per cell in horizontal direction.
- `nmarkers_cell_y`: Number of markers per cell in vertical direction.
- `mynum`: Total number of markers in vertical direction.
- `mxnum`: Total number of markers in horizontal direction.
- `marknum`: Total number of markers.
- `mystep`: Distance between markers in vertical direction in meters.
- `mxstep`: Distance between markers in horizontal direction in meters.

# Returns
- A new Distribution struct with the given marker parameters.
"""
mutable struct Distribution <: AbstractParameterGroup
    iuse_random::ParameterInt
    dx_marker::ParameterFloat
    dy_marker::ParameterFloat
    nmarkers_cell_x::ParameterFloat
    nmarkers_cell_y::ParameterFloat
    mxnum::ParameterInt
    mynum::ParameterInt
    marknum::ParameterInt
    mxstep::ParameterFloat
    mystep::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Distribution(
    marker_parameters::NamedTuple
)::Distribution
    dx_marker = marker_parameters.dx_marker
    dy_marker = marker_parameters.dy_marker
    nmarkers_cell_y = marker_parameters.nmarkers_cell_y
    nmarkers_cell_x = marker_parameters.nmarkers_cell_x
    mynum = marker_parameters.mynum
    mxnum = marker_parameters.mxnum
    marknum = marker_parameters.marknum
    mystep = marker_parameters.mystep
    mxstep = marker_parameters.mxstep

    pdata = get_eb_parameters()
    data = Distribution(
        ParameterInt(0, pdata.iuse_random.name, pdata.iuse_random.units, pdata.iuse_random.description),
        ParameterFloat(dx_marker, pdata.dx_marker.name, pdata.dx_marker.units, pdata.dx_marker.description),
        ParameterFloat(dy_marker, pdata.dy_marker.name, pdata.dy_marker.units, pdata.dy_marker.description),
        ParameterFloat(nmarkers_cell_x, pdata.nmarkers_cell_x.name, pdata.nmarkers_cell_x.units, pdata.nmarkers_cell_x.description),
        ParameterFloat(nmarkers_cell_y, pdata.nmarkers_cell_y.name, pdata.nmarkers_cell_y.units, pdata.nmarkers_cell_y.description),
        ParameterInt(mxnum, pdata.mxnum.name, pdata.mxnum.units, pdata.mxnum.description),
        ParameterInt(mynum, pdata.mynum.name, pdata.mynum.units, pdata.mynum.description),
        ParameterInt(marknum, pdata.marknum.name, pdata.marknum.units, pdata.marknum.description),
        ParameterFloat(mxstep, pdata.mxstep.name, pdata.mxstep.units, pdata.mxstep.description),
        ParameterFloat(mystep, pdata.mystep.name, pdata.mystep.units, pdata.mystep.description),
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
