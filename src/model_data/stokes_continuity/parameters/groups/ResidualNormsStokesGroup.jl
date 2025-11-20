module ResidualNormsStokesGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.stokes_continuity.parameters"
const GRP_NAME = "residual_norms"

const PDATA = get_eb_parameters()

"""
    ResidualNormsStokes <: AbstractParameterGroup

Parameter group for Stokes residual norm tracking.

# Fields
- `resnl_L2_ini::`[`ParameterFloat`](@ref): $(PDATA.resnl_L2_ini.description)
- `resnl_L2::`[`ParameterFloat`](@ref): $(PDATA.resnl_L2.description)
- `resnl_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.resnl_rel_L2.description)
- `resx_L2::`[`ParameterFloat`](@ref): $(PDATA.resx_L2.description)
- `resy_L2::`[`ParameterFloat`](@ref): $(PDATA.resy_L2.description)
- `resc_L2::`[`ParameterFloat`](@ref): $(PDATA.resc_L2.description)
- `resx_L2_ini::`[`ParameterFloat`](@ref): $(PDATA.resx_L2_ini.description)
- `resy_L2_ini::`[`ParameterFloat`](@ref): $(PDATA.resy_L2_ini.description)
- `resc_L2_ini::`[`ParameterFloat`](@ref): $(PDATA.resc_L2_ini.description)
- `resnlx_L2::`[`ParameterFloat`](@ref): $(PDATA.resnlx_L2.description)
- `resnly_L2::`[`ParameterFloat`](@ref): $(PDATA.resnly_L2.description)
- `resnlc_L2::`[`ParameterFloat`](@ref): $(PDATA.resnlc_L2.description)
- `resnlx_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.resnlx_rel_L2.description)
- `resnly_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.resnly_rel_L2.description)
- `resnlc_rel_L2::`[`ParameterFloat`](@ref): $(PDATA.resnlc_rel_L2.description)
- `resnlx_L2_ini::`[`ParameterFloat`](@ref): $(PDATA.resnlx_L2_ini.description)
- `resnly_L2_ini::`[`ParameterFloat`](@ref): $(PDATA.resnly_L2_ini.description)
- `resnlc_L2_ini::`[`ParameterFloat`](@ref): $(PDATA.resnlc_L2_ini.description)
- `obj_list::Vector{ParameterFloat}`: List of parameter objects

# Nested Dot Access
- `resnl_L2_ini = $(ROOT_NAME).$(GRP_NAME).resnl_L2_ini.value`
- `resnl_L2 = $(ROOT_NAME).$(GRP_NAME).resnl_L2.value`
- `resnl_rel_L2 = $(ROOT_NAME).$(GRP_NAME).resnl_rel_L2.value`
- `resx_L2 = $(ROOT_NAME).$(GRP_NAME).resx_L2.value`
- `resy_L2 = $(ROOT_NAME).$(GRP_NAME).resy_L2.value`
- `resc_L2 = $(ROOT_NAME).$(GRP_NAME).resc_L2.value`

# Constructor
    ResidualNormsStokes()

# Returns
- `ResidualNormsStokes`: New ResidualNormsStokes parameter group with initialized values

"""
mutable struct ResidualNormsStokes <: AbstractParameterGroup
    resnl_L2_ini::ParameterFloat
    resnl_L2::ParameterFloat
    resnl_rel_L2::ParameterFloat
    resx_L2::ParameterFloat
    resy_L2::ParameterFloat
    resc_L2::ParameterFloat
    resx_L2_ini::ParameterFloat
    resy_L2_ini::ParameterFloat
    resc_L2_ini::ParameterFloat
    resnlx_L2::ParameterFloat
    resnly_L2::ParameterFloat
    resnlc_L2::ParameterFloat
    resnlx_rel_L2::ParameterFloat
    resnly_rel_L2::ParameterFloat
    resnlc_rel_L2::ParameterFloat
    resnlx_L2_ini::ParameterFloat
    resnly_L2_ini::ParameterFloat
    resnlc_L2_ini::ParameterFloat
    obj_list::Vector{ParameterFloat}
end

function ResidualNormsStokes()::ResidualNormsStokes
    pdata = get_eb_parameters()
    data = ResidualNormsStokes(
        pdata.resnl_L2_ini,
        pdata.resnl_L2,
        pdata.resnl_rel_L2,
        pdata.resx_L2,
        pdata.resy_L2,
        pdata.resc_L2,
        pdata.resx_L2_ini,
        pdata.resy_L2_ini,
        pdata.resc_L2_ini,
        pdata.resnlx_L2,
        pdata.resnly_L2,
        pdata.resnlc_L2,
        pdata.resnlx_rel_L2,
        pdata.resnly_rel_L2,
        pdata.resnlc_rel_L2,
        pdata.resnlx_L2_ini,
        pdata.resnly_L2_ini,
        pdata.resnlc_L2_ini,
        ParameterFloat[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
