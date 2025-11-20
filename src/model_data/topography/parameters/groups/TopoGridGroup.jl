module TopoGridGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.topography.parameters"
const GRP_NAME = "topo_grid"

const PDATA = get_eb_parameters()

"""
    TopoGrid <: AbstractParameterGroup

Parameter group for topography grid parameters.

# Fields
- `topo_xsize::`[`ParameterFloat`](@ref): $(PDATA.topo_xsize.description)
- `toponum::`[`ParameterInt`](@ref): $(PDATA.toponum.description)
- `dx_topo::`[`ParameterFloat`](@ref): $(PDATA.dx_topo.description)
- `nsmooth_top_bottom::`[`ParameterInt`](@ref): $(PDATA.nsmooth_top_bottom.description)
- `marker_search_factor::`[`ParameterFloat`](@ref): $(PDATA.marker_search_factor.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `topo_xsize = $(ROOT_NAME).$(GRP_NAME).topo_xsize.value`
- `toponum = $(ROOT_NAME).$(GRP_NAME).toponum.value`
- `dx_topo = $(ROOT_NAME).$(GRP_NAME).dx_topo.value`
- `nsmooth_top_bottom = $(ROOT_NAME).$(GRP_NAME).nsmooth_top_bottom.value`
- `marker_search_factor = $(ROOT_NAME).$(GRP_NAME).marker_search_factor.value`

# Constructor
    TopoGrid()

# Returns
- `TopoGrid`: New TopoGrid parameter group with initialized values

"""
mutable struct TopoGrid <: AbstractParameterGroup
    topo_xsize::ParameterFloat
    toponum::ParameterInt
    dx_topo::ParameterFloat
    nsmooth_top_bottom::ParameterInt
    marker_search_factor::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function TopoGrid()::TopoGrid
    pdata = get_eb_parameters()
    data = TopoGrid(
        pdata.topo_xsize,
        pdata.toponum,
        pdata.dx_topo,
        pdata.nsmooth_top_bottom,
        pdata.marker_search_factor,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
