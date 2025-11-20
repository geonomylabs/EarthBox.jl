module GridGeometryGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.grids.parameters"
const GRP_NAME = "geometry"

const PDATA = get_eb_parameters()

"""
    GridGeometry <: AbstractParameterGroup

Parameter group for basic grid geometry parameters.

# Fields
- `ynum::`[`ParameterInt`](@ref): $(PDATA.ynum.description)
- `xnum::`[`ParameterInt`](@ref): $(PDATA.xnum.description)
- `ysize::`[`ParameterFloat`](@ref): $(PDATA.ysize.description)
- `xsize::`[`ParameterFloat`](@ref): $(PDATA.xsize.description)
- `ystpavg::`[`ParameterFloat`](@ref): $(PDATA.ystpavg.description)
- `xstpavg::`[`ParameterFloat`](@ref): $(PDATA.xstpavg.description)
- `ymin::`[`ParameterFloat`](@ref): $(PDATA.ymin.description)
- `ymax::`[`ParameterFloat`](@ref): $(PDATA.ymax.description)
- `xmin::`[`ParameterFloat`](@ref): $(PDATA.xmin.description)
- `xmax::`[`ParameterFloat`](@ref): $(PDATA.xmax.description)
- `xsize_start::`[`ParameterFloat`](@ref): $(PDATA.xsize_start.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `ynum = $(ROOT_NAME).$(GRP_NAME).ynum.value`
- `xnum = $(ROOT_NAME).$(GRP_NAME).xnum.value`
- `ysize = $(ROOT_NAME).$(GRP_NAME).ysize.value`
- `xsize = $(ROOT_NAME).$(GRP_NAME).xsize.value`
- `ystpavg = $(ROOT_NAME).$(GRP_NAME).ystpavg.value`
- `xstpavg = $(ROOT_NAME).$(GRP_NAME).xstpavg.value`
- `ymin = $(ROOT_NAME).$(GRP_NAME).ymin.value`
- `ymax = $(ROOT_NAME).$(GRP_NAME).ymax.value`
- `xmin = $(ROOT_NAME).$(GRP_NAME).xmin.value`
- `xmax = $(ROOT_NAME).$(GRP_NAME).xmax.value`
- `xsize_start = $(ROOT_NAME).$(GRP_NAME).xsize_start.value`

# Constructor
    GridGeometry(ynum::Int, xnum::Int, ysize::Float64, xsize::Float64)

Create a new GridGeometry parameter group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction
- `ysize::Float64`: Domain size in y-direction in meters
- `xsize::Float64`: Domain size in x-direction in meters

# Returns
- `GridGeometry`: New GridGeometry parameter group with initialized values

"""
mutable struct GridGeometry <: AbstractParameterGroup
    ynum::ParameterInt
    xnum::ParameterInt
    ysize::ParameterFloat
    xsize::ParameterFloat
    ystpavg::ParameterFloat
    xstpavg::ParameterFloat
    ymin::ParameterFloat
    ymax::ParameterFloat
    xmin::ParameterFloat
    xmax::ParameterFloat
    xsize_start::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function GridGeometry(
    ynum::Int, xnum::Int, ysize::Float64, xsize::Float64
)::GridGeometry
    pdata = get_eb_parameters()
    data = GridGeometry(
        ParameterInt(ynum, pdata.ynum.name, pdata.ynum.units, pdata.ynum.description),
        ParameterInt(xnum, pdata.xnum.name, pdata.xnum.units, pdata.xnum.description),
        ParameterFloat(ysize, pdata.ysize.name, pdata.ysize.units, pdata.ysize.description),
        ParameterFloat(xsize, pdata.xsize.name, pdata.xsize.units, pdata.xsize.description),
        pdata.ystpavg,
        pdata.xstpavg,
        pdata.ymin,
        ParameterFloat(ysize, pdata.ymax.name, pdata.ymax.units, pdata.ymax.description),
        pdata.xmin,
        ParameterFloat(xsize, pdata.xmax.name, pdata.xmax.units, pdata.xmax.description),
        ParameterFloat(xsize, pdata.xsize_start.name, pdata.xsize_start.units, pdata.xsize_start.description),
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 