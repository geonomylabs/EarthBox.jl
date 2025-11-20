module BuildSysStokesGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.stokes_continuity.parameters"
const GRP_NAME = "build"

const PDATA = get_eb_parameters()

"""
    BuildSysStokes <: AbstractParameterGroup

Parameter group for sparse matrix build parameters for Stokes system.

# Fields
- `ibuild::`[`ParameterInt`](@ref): $(PDATA.ibuild_stokes.description)
- `hshift_to_vxR::`[`ParameterInt`](@ref): $(PDATA.hshift_to_vxR.description)
- `N::`[`ParameterInt`](@ref): $(PDATA.Nstokes.description)
- `nonzero_max::`[`ParameterInt`](@ref): $(PDATA.nonzero_max_stokes.description)
- `pscale::`[`ParameterFloat`](@ref): $(PDATA.pscale.description)
- `iuse_interface_stabilization::`[`ParameterInt`](@ref): $(PDATA.iuse_interface_stabilization.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `ibuild = $(ROOT_NAME).$(GRP_NAME).ibuild.value`
- `hshift_to_vxR = $(ROOT_NAME).$(GRP_NAME).hshift_to_vxR.value`
- `N = $(ROOT_NAME).$(GRP_NAME).N.value`
- `nonzero_max = $(ROOT_NAME).$(GRP_NAME).nonzero_max.value`
- `pscale = $(ROOT_NAME).$(GRP_NAME).pscale.value`
- `iuse_interface_stabilization = $(ROOT_NAME).$(GRP_NAME).iuse_interface_stabilization.value`

# Constructor
    BuildSysStokes(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `BuildSysStokes`: New BuildSysStokes parameter group with initialized values

"""
mutable struct BuildSysStokes <: AbstractParameterGroup
    ibuild::ParameterInt
    hshift_to_vxR::ParameterInt
    N::ParameterInt
    nonzero_max::ParameterInt
    pscale::ParameterFloat
    iuse_interface_stabilization::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function BuildSysStokes(ynum::Int, xnum::Int)::BuildSysStokes
    pdata = get_eb_parameters()
    data = BuildSysStokes(
        ParameterInt(1, pdata.ibuild_stokes.name, pdata.ibuild_stokes.units, pdata.ibuild_stokes.description),
        ParameterInt((ynum-1)*3, pdata.hshift_to_vxR.name, pdata.hshift_to_vxR.units, pdata.hshift_to_vxR.description),
        ParameterInt((xnum-1)*(ynum-1)*3, pdata.Nstokes.name, pdata.Nstokes.units, pdata.Nstokes.description),
        ParameterInt(xnum*ynum*31, pdata.nonzero_max_stokes.name, pdata.nonzero_max_stokes.units, pdata.nonzero_max_stokes.description),
        ParameterFloat(1.0, pdata.pscale.name, pdata.pscale.units, pdata.pscale.description),
        ParameterInt(0, pdata.iuse_interface_stabilization.name, pdata.iuse_interface_stabilization.units, pdata.iuse_interface_stabilization.description),
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
