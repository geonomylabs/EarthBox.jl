"""
Module for velocity boundary condition parameters.

Provides data structures for configuring velocity boundary conditions including
extension, contraction, shear, rotation, and inflow velocities.
"""
module VelocityGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: AbstractParameterGroup
import EarthBox.Parameters: ParameterFloat, ParameterInt

const ROOT_NAME = "model.bcs.parameters"
const GRP_NAME = "velocity"

const PDATA = get_eb_parameters()

"""
    Velocity <: AbstractParameterGroup

Parameter group for velocity boundary conditions.

# Fields
- `velocity::`[`ParameterFloat`](@ref): $(PDATA.velocity.description)
- `full_velocity_extension::`[`ParameterFloat`](@ref): $(PDATA.full_velocity_extension.description)
- `full_velocity_extension_step1::`[`ParameterFloat`](@ref): $(PDATA.full_velocity_extension_step1.description)
- `full_velocity_extension_step2::`[`ParameterFloat`](@ref): $(PDATA.full_velocity_extension_step2.description)
- `full_velocity_contraction::`[`ParameterFloat`](@ref): $(PDATA.full_velocity_contraction.description)
- `velocity_shear::`[`ParameterFloat`](@ref): $(PDATA.velocity_shear.description)
- `velocity_rotation::`[`ParameterFloat`](@ref): $(PDATA.velocity_rotation.description)
- `iuse_strain_rate::`[`ParameterInt`](@ref): $(PDATA.iuse_strain_rate.description)
- `strain_rate_bc::`[`ParameterFloat`](@ref): $(PDATA.strain_rate_bc.description)
- `vyu::`[`ParameterFloat`](@ref): $(PDATA.vyu.description)
- `vyl::`[`ParameterFloat`](@ref): $(PDATA.vyl.description)
- `velocity_internal_x::`[`ParameterFloat`](@ref): $(PDATA.velocity_internal_x.description)
- `velocity_internal_y::`[`ParameterFloat`](@ref): $(PDATA.velocity_internal_y.description)
- `plate_thickness::`[`ParameterFloat`](@ref): $(PDATA.plate_thickness.description)
- `smoothing_thickness::`[`ParameterFloat`](@ref): $(PDATA.smoothing_thickness.description)
- `velocity_inflow_left::`[`ParameterFloat`](@ref): $(PDATA.velocity_inflow_left.description)
- `velocity_inflow_right::`[`ParameterFloat`](@ref): $(PDATA.velocity_inflow_right.description)
- `velocity_inflow_smooth_avg_left::`[`ParameterFloat`](@ref): $(PDATA.velocity_inflow_smooth_avg_left.description)
- `velocity_inflow_smooth_avg_right::`[`ParameterFloat`](@ref): $(PDATA.velocity_inflow_smooth_avg_right.description)
- `iuse_velocity_stop::`[`ParameterInt`](@ref): $(PDATA.iuse_velocity_stop.description)
- `velocity_stop_time::`[`ParameterFloat`](@ref): $(PDATA.velocity_stop_time.description)
- `ivelocity_stop_counter::`[`ParameterInt`](@ref): $(PDATA.ivelocity_stop_counter.description)

# Nested Dot Access
- `velocity = $(ROOT_NAME).$(GRP_NAME).velocity.value`
- `full_velocity_extension = $(ROOT_NAME).$(GRP_NAME).full_velocity_extension.value`
- `full_velocity_extension_step1 = $(ROOT_NAME).$(GRP_NAME).full_velocity_extension_step1.value`
- `full_velocity_extension_step2 = $(ROOT_NAME).$(GRP_NAME).full_velocity_extension_step2.value`
- `full_velocity_contraction = $(ROOT_NAME).$(GRP_NAME).full_velocity_contraction.value`
- `velocity_shear = $(ROOT_NAME).$(GRP_NAME).velocity_shear.value`
- `velocity_rotation = $(ROOT_NAME).$(GRP_NAME).velocity_rotation.value`
- `iuse_strain_rate = $(ROOT_NAME).$(GRP_NAME).iuse_strain_rate.value`
- `strain_rate_bc = $(ROOT_NAME).$(GRP_NAME).strain_rate_bc.value`
- `vyu = $(ROOT_NAME).$(GRP_NAME).vyu.value`
- `vyl = $(ROOT_NAME).$(GRP_NAME).vyl.value`
- `velocity_internal_x = $(ROOT_NAME).$(GRP_NAME).velocity_internal_x.value`
- `velocity_internal_y = $(ROOT_NAME).$(GRP_NAME).velocity_internal_y.value`
- `plate_thickness = $(ROOT_NAME).$(GRP_NAME).plate_thickness.value`
- `smoothing_thickness = $(ROOT_NAME).$(GRP_NAME).smoothing_thickness.value`
- `velocity_inflow_left = $(ROOT_NAME).$(GRP_NAME).velocity_inflow_left.value`
- `velocity_inflow_right = $(ROOT_NAME).$(GRP_NAME).velocity_inflow_right.value`
- `velocity_inflow_smooth_avg_left = $(ROOT_NAME).$(GRP_NAME).velocity_inflow_smooth_avg_left.value`
- `velocity_inflow_smooth_avg_right = $(ROOT_NAME).$(GRP_NAME).velocity_inflow_smooth_avg_right.value`
- `iuse_velocity_stop = $(ROOT_NAME).$(GRP_NAME).iuse_velocity_stop.value`
- `velocity_stop_time = $(ROOT_NAME).$(GRP_NAME).velocity_stop_time.value`
- `ivelocity_stop_counter = $(ROOT_NAME).$(GRP_NAME).ivelocity_stop_counter.value`

# Constructor
    Velocity()

Initializes velocity boundary condition parameters with default values.
"""
mutable struct Velocity <: AbstractParameterGroup
    velocity::ParameterFloat
    full_velocity_extension::ParameterFloat
    full_velocity_extension_step1::ParameterFloat
    full_velocity_extension_step2::ParameterFloat
    full_velocity_contraction::ParameterFloat
    velocity_shear::ParameterFloat
    velocity_rotation::ParameterFloat
    iuse_strain_rate::ParameterInt
    strain_rate_bc::ParameterFloat
    vyu::ParameterFloat
    vyl::ParameterFloat
    velocity_internal_x::ParameterFloat
    velocity_internal_y::ParameterFloat
    plate_thickness::ParameterFloat
    smoothing_thickness::ParameterFloat
    velocity_inflow_left::ParameterFloat
    velocity_inflow_right::ParameterFloat
    velocity_inflow_smooth_avg_left::ParameterFloat
    velocity_inflow_smooth_avg_right::ParameterFloat
    iuse_velocity_stop::ParameterInt
    velocity_stop_time::ParameterFloat
    ivelocity_stop_counter::ParameterInt
end

function Velocity()::Velocity
    pdata = get_eb_parameters()
    data = Velocity(
        pdata.velocity,
        pdata.full_velocity_extension,
        pdata.full_velocity_extension_step1,
        pdata.full_velocity_extension_step2,
        pdata.full_velocity_contraction,
        pdata.velocity_shear,
        pdata.velocity_rotation,
        pdata.iuse_strain_rate,
        pdata.strain_rate_bc,
        pdata.vyu,
        pdata.vyl,
        pdata.velocity_internal_x,
        pdata.velocity_internal_y,
        pdata.plate_thickness,
        pdata.smoothing_thickness,
        pdata.velocity_inflow_left,
        pdata.velocity_inflow_right,
        pdata.velocity_inflow_smooth_avg_left,
        pdata.velocity_inflow_smooth_avg_right,
        pdata.iuse_velocity_stop,
        pdata.velocity_stop_time,
        pdata.ivelocity_stop_counter
    )
    return data
end

end # module
