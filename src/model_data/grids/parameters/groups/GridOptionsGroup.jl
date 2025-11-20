module GridOptionsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterInt, ParameterStr
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.grids.parameters"
const GRP_NAME = "grid_options"

const PDATA = get_eb_parameters()

"""
    GridOptions <: AbstractParameterGroup

Parameter group for grid options.

# Fields
- `itype_grid::`[`ParameterInt`](@ref): $(PDATA.itype_grid.description)
- `stype_grid::`[`ParameterStr`](@ref): $(PDATA.stype_grid.description)
- `obj_list::Vector{ParameterInt}`: List of parameter objects

# Nested Dot Access
- `itype_grid = $(ROOT_NAME).$(GRP_NAME).itype_grid.value`
- `stype_grid = $(ROOT_NAME).$(GRP_NAME).stype_grid.value`

# Constructor
    GridOptions()

Create a new GridOptions parameter group with default values.

# Returns
- `GridOptions`: New GridOptions parameter group with initialized values

"""
mutable struct GridOptions <: AbstractParameterGroup
    itype_grid::ParameterInt
    stype_grid::ParameterStr
    obj_list::Vector{ParameterInt}
end

function GridOptions()::GridOptions
    pdata = get_eb_parameters()
    data = GridOptions(
        pdata.itype_grid,
        pdata.stype_grid,
        Union{ParameterInt, ParameterStr}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module 