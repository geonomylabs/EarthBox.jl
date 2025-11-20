module BuildGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.heat_equation.parameters"
const GRP_NAME = "build"

const PDATA = get_eb_parameters()

"""
    Build <: AbstractParameterGroup

Parameter group for sparse matrix build parameters.

# Fields
- `ibuild::`[`ParameterInt`](@ref): $(PDATA.ibuild_heat.description)
- `N::`[`ParameterInt`](@ref): $(PDATA.Nheat.description)
- `nonzero_max_heat::`[`ParameterInt`](@ref): $(PDATA.nonzero_max_heat.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `ibuild = $(ROOT_NAME).$(GRP_NAME).ibuild.value`
- `N = $(ROOT_NAME).$(GRP_NAME).N.value`
- `nonzero_max_heat = $(ROOT_NAME).$(GRP_NAME).nonzero_max_heat.value`

# Constructor
    Build(xnum::Int, ynum::Int)

# Arguments
- `xnum::Int`: Number of grid points in x-direction
- `ynum::Int`: Number of grid points in y-direction

# Returns
- `Build`: New Build parameter group with initialized values

"""
mutable struct Build <: AbstractParameterGroup
    ibuild::ParameterInt
    N::ParameterInt
    nonzero_max_heat::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Build(xnum::Int, ynum::Int)::Build
    pdata = get_eb_parameters()
    data = Build(
        ParameterInt(1, pdata.ibuild_heat.name, pdata.ibuild_heat.units, pdata.ibuild_heat.description),
        ParameterInt(xnum*ynum, pdata.Nheat.name, pdata.Nheat.units, pdata.Nheat.description),
        ParameterInt(xnum*ynum*31, pdata.nonzero_max_heat.name, pdata.nonzero_max_heat.units, pdata.nonzero_max_heat.description),
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
