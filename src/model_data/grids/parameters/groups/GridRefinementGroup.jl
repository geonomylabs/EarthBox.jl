module GridRefinementGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.grids.parameters"
const GRP_NAME = "refinement"

const PDATA = get_eb_parameters()

"""
    GridRefinement <: AbstractParameterGroup

Parameter group for grid refinement parameters.

# Fields
- `dx_highres::`[`ParameterFloat`](@ref): $(PDATA.dx_highres.description)
- `dx_lowres::`[`ParameterFloat`](@ref): $(PDATA.dx_lowres.description)
- `xo_highres::`[`ParameterFloat`](@ref): $(PDATA.xo_highres.description)
- `ixo_highres::`[`ParameterInt`](@ref): $(PDATA.ixo_highres.description)
- `xf_highres::`[`ParameterFloat`](@ref): $(PDATA.xf_highres.description)
- `dy_highres::`[`ParameterFloat`](@ref): $(PDATA.dy_highres.description)
- `dy_lowres::`[`ParameterFloat`](@ref): $(PDATA.dy_lowres.description)
- `yf_highres::`[`ParameterFloat`](@ref): $(PDATA.yf_highres.description)
- `iuse_trench::`[`ParameterInt`](@ref): $(PDATA.iuse_trench.description)
- `trench_location::`[`ParameterFloat`](@ref): $(PDATA.trench_location.description)
- `iuse_refinement_delay::`[`ParameterInt`](@ref): $(PDATA.iuse_refinement_delay.description)
- `refinement_time::`[`ParameterFloat`](@ref): $(PDATA.refinement_time.description)
- `refinement_flag::`[`ParameterInt`](@ref): $(PDATA.refinement_flag.description)
- `iuse_refinement_gap::`[`ParameterInt`](@ref): $(PDATA.iuse_refinement_gap.description)
- `refinement_gap_start_time::`[`ParameterFloat`](@ref): $(PDATA.refinement_gap_start_time.description)
- `refinement_gap_end_time::`[`ParameterFloat`](@ref): $(PDATA.refinement_gap_end_time.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `dx_highres = $(ROOT_NAME).$(GRP_NAME).dx_highres.value`
- `dx_lowres = $(ROOT_NAME).$(GRP_NAME).dx_lowres.value`
- `xo_highres = $(ROOT_NAME).$(GRP_NAME).xo_highres.value`
- `ixo_highres = $(ROOT_NAME).$(GRP_NAME).ixo_highres.value`
- `xf_highres = $(ROOT_NAME).$(GRP_NAME).xf_highres.value`
- `dy_highres = $(ROOT_NAME).$(GRP_NAME).dy_highres.value`
- `dy_lowres = $(ROOT_NAME).$(GRP_NAME).dy_lowres.value`
- `yf_highres = $(ROOT_NAME).$(GRP_NAME).yf_highres.value`
- `iuse_trench = $(ROOT_NAME).$(GRP_NAME).iuse_trench.value`
- `trench_location = $(ROOT_NAME).$(GRP_NAME).trench_location.value`
- `iuse_refinement_delay = $(ROOT_NAME).$(GRP_NAME).iuse_refinement_delay.value`
- `refinement_time = $(ROOT_NAME).$(GRP_NAME).refinement_time.value`
- `refinement_flag = $(ROOT_NAME).$(GRP_NAME).refinement_flag.value`
- `iuse_refinement_gap = $(ROOT_NAME).$(GRP_NAME).iuse_refinement_gap.value`
- `refinement_gap_start_time = $(ROOT_NAME).$(GRP_NAME).refinement_gap_start_time.value`
- `refinement_gap_end_time = $(ROOT_NAME).$(GRP_NAME).refinement_gap_end_time.value`

# Constructor
    GridRefinement()

Create a new GridRefinement parameter group with default values.

# Returns
- `GridRefinement`: New GridRefinement parameter group with initialized values

"""
mutable struct GridRefinement <: AbstractParameterGroup
    dx_highres::ParameterFloat
    dx_lowres::ParameterFloat
    xo_highres::ParameterFloat
    ixo_highres::ParameterInt
    xf_highres::ParameterFloat
    dy_highres::ParameterFloat
    dy_lowres::ParameterFloat
    yf_highres::ParameterFloat
    iuse_trench::ParameterInt
    trench_location::ParameterFloat
    iuse_refinement_delay::ParameterInt
    refinement_time::ParameterFloat
    refinement_flag::ParameterInt
    iuse_refinement_gap::ParameterInt
    refinement_gap_start_time::ParameterFloat
    refinement_gap_end_time::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function GridRefinement()::GridRefinement
    pdata = get_eb_parameters()
    # Negative values NaN and -9999 indicate that parameter has not been set.
    data = GridRefinement(
        pdata.dx_highres,
        pdata.dx_lowres,
        pdata.xo_highres,
        pdata.ixo_highres,
        pdata.xf_highres,
        pdata.dy_highres,
        pdata.dy_lowres,
        pdata.yf_highres,
        pdata.iuse_trench,
        pdata.trench_location,
        pdata.iuse_refinement_delay,
        pdata.refinement_time,
        pdata.refinement_flag,
        pdata.iuse_refinement_gap,
        pdata.refinement_gap_start_time,
        pdata.refinement_gap_end_time,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 