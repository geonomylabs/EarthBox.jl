module ExtrusionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.melting.parameters"
const GRP_NAME = "extrusion"

const PDATA = get_eb_parameters()

"""
    Extrusion <: AbstractParameterGroup

Parameter group for melt extrusion model parameters.

# Fields
- `iuse_extrusion::`[`ParameterInt`](@ref): $(PDATA.iuse_extrusion.description)
- `extrusion_volume_factor::`[`ParameterFloat`](@ref): $(PDATA.extrusion_volume_factor.description)
- `extrusion_volume_factor_max::`[`ParameterFloat`](@ref): $(PDATA.extrusion_volume_factor_max.description)
- `characteristic_magmatic_crust_height::`[`ParameterFloat`](@ref): $(PDATA.characteristic_magmatic_crust_height.description)
- `characteristic_magmatic_crust_height_min::`[`ParameterFloat`](@ref): $(PDATA.characteristic_magmatic_crust_height_min.description)
- `characteristic_magmatic_crust_height_max::`[`ParameterFloat`](@ref): $(PDATA.characteristic_magmatic_crust_height_max.description)
- `characteristic_flow_length_subaerial::`[`ParameterFloat`](@ref): $(PDATA.characteristic_flow_length_subaerial.description)
- `characteristic_flow_length_submarine::`[`ParameterFloat`](@ref): $(PDATA.characteristic_flow_length_submarine.description)
- `residual_lava_thickness_subaerial::`[`ParameterFloat`](@ref): $(PDATA.residual_lava_thickness_subaerial.description)
- `residual_lava_thickness_submarine::`[`ParameterFloat`](@ref): $(PDATA.residual_lava_thickness_submarine.description)
- `iuse_random_eruption_location::`[`ParameterInt`](@ref): $(PDATA.iuse_random_eruption_location.description)
- `iuse_normal_eruption_location::`[`ParameterInt`](@ref): $(PDATA.iuse_normal_eruption_location.description)
- `decimation_factor::`[`ParameterInt`](@ref): $(PDATA.decimation_factor.description)
- `initial_magma_flush_steps::`[`ParameterInt`](@ref): $(PDATA.initial_magma_flush_steps.description)
- `magma_flush_factor::`[`ParameterFloat`](@ref): $(PDATA.magma_flush_factor.description)
- `width_eruption_domain_fixed::`[`ParameterFloat`](@ref): $(PDATA.width_eruption_domain_fixed.description)
- `width_eruption_domain_fixed_max::`[`ParameterFloat`](@ref): $(PDATA.width_eruption_domain_fixed_max.description)
- `porosity_initial_lava_flow::`[`ParameterFloat`](@ref): $(PDATA.porosity_initial_lava_flow.description)
- `decay_depth_lava_flow::`[`ParameterFloat`](@ref): $(PDATA.decay_depth_lava_flow.description)
- `time_of_next_eruption_myr::`[`ParameterFloat`](@ref): $(PDATA.time_of_next_eruption_myr.description)
- `eruption_interval_yr::`[`ParameterFloat`](@ref): $(PDATA.eruption_interval_yr.description)
- `iuse_eruption_interval::`[`ParameterInt`](@ref): $(PDATA.iuse_eruption_interval.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `iuse_extrusion = $(ROOT_NAME).$(GRP_NAME).iuse_extrusion.value`
- `extrusion_volume_factor = $(ROOT_NAME).$(GRP_NAME).extrusion_volume_factor.value`
- `extrusion_volume_factor_max = $(ROOT_NAME).$(GRP_NAME).extrusion_volume_factor_max.value`
- `characteristic_magmatic_crust_height = $(ROOT_NAME).$(GRP_NAME).characteristic_magmatic_crust_height.value`
- `characteristic_magmatic_crust_height_min = $(ROOT_NAME).$(GRP_NAME).characteristic_magmatic_crust_height_min.value`
- `characteristic_magmatic_crust_height_max = $(ROOT_NAME).$(GRP_NAME).characteristic_magmatic_crust_height_max.value`
- `characteristic_flow_length_subaerial = $(ROOT_NAME).$(GRP_NAME).characteristic_flow_length_subaerial.value`
- `characteristic_flow_length_submarine = $(ROOT_NAME).$(GRP_NAME).characteristic_flow_length_submarine.value`
- `residual_lava_thickness_subaerial = $(ROOT_NAME).$(GRP_NAME).residual_lava_thickness_subaerial.value`
- `residual_lava_thickness_submarine = $(ROOT_NAME).$(GRP_NAME).residual_lava_thickness_submarine.value`
- `iuse_random_eruption_location = $(ROOT_NAME).$(GRP_NAME).iuse_random_eruption_location.value`
- `iuse_normal_eruption_location = $(ROOT_NAME).$(GRP_NAME).iuse_normal_eruption_location.value`
- `decimation_factor = $(ROOT_NAME).$(GRP_NAME).decimation_factor.value`
- `initial_magma_flush_steps = $(ROOT_NAME).$(GRP_NAME).initial_magma_flush_steps.value`
- `magma_flush_factor = $(ROOT_NAME).$(GRP_NAME).magma_flush_factor.value`
- `width_eruption_domain_fixed = $(ROOT_NAME).$(GRP_NAME).width_eruption_domain_fixed.value`
- `width_eruption_domain_fixed_max = $(ROOT_NAME).$(GRP_NAME).width_eruption_domain_fixed_max.value`
- `porosity_initial_lava_flow = $(ROOT_NAME).$(GRP_NAME).porosity_initial_lava_flow.value`
- `decay_depth_lava_flow = $(ROOT_NAME).$(GRP_NAME).decay_depth_lava_flow.value`
- `time_of_next_eruption_myr = $(ROOT_NAME).$(GRP_NAME).time_of_next_eruption_myr.value`
- `eruption_interval_yr = $(ROOT_NAME).$(GRP_NAME).eruption_interval_yr.value`
- `iuse_eruption_interval = $(ROOT_NAME).$(GRP_NAME).iuse_eruption_interval.value`

# Constructor
    Extrusion()

Create a new Extrusion parameter group with default values.

# Returns
- `Extrusion`: New Extrusion parameter group with initialized values

"""
mutable struct Extrusion <: AbstractParameterGroup
    iuse_extrusion::ParameterInt
    extrusion_volume_factor::ParameterFloat
    extrusion_volume_factor_max::ParameterFloat
    characteristic_magmatic_crust_height::ParameterFloat
    characteristic_magmatic_crust_height_min::ParameterFloat
    characteristic_magmatic_crust_height_max::ParameterFloat
    characteristic_flow_length_subaerial::ParameterFloat
    characteristic_flow_length_submarine::ParameterFloat
    residual_lava_thickness_subaerial::ParameterFloat
    residual_lava_thickness_submarine::ParameterFloat
    iuse_random_eruption_location::ParameterInt
    iuse_normal_eruption_location::ParameterInt
    decimation_factor::ParameterInt
    initial_magma_flush_steps::ParameterInt
    magma_flush_factor::ParameterFloat
    width_eruption_domain_fixed::ParameterFloat
    width_eruption_domain_fixed_max::ParameterFloat
    porosity_initial_lava_flow::ParameterFloat
    decay_depth_lava_flow::ParameterFloat
    time_of_next_eruption_myr::ParameterFloat
    eruption_interval_yr::ParameterFloat
    iuse_eruption_interval::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Extrusion()::Extrusion
    pdata = get_eb_parameters()
    data = Extrusion(
        pdata.iuse_extrusion,
        pdata.extrusion_volume_factor,
        pdata.extrusion_volume_factor_max,
        pdata.characteristic_magmatic_crust_height,
        pdata.characteristic_magmatic_crust_height_min,
        pdata.characteristic_magmatic_crust_height_max,
        pdata.characteristic_flow_length_subaerial,
        pdata.characteristic_flow_length_submarine,
        pdata.residual_lava_thickness_subaerial,
        pdata.residual_lava_thickness_submarine,
        pdata.iuse_random_eruption_location,
        pdata.iuse_normal_eruption_location,
        pdata.decimation_factor,
        pdata.initial_magma_flush_steps,
        pdata.magma_flush_factor,
        pdata.width_eruption_domain_fixed,
        pdata.width_eruption_domain_fixed_max,
        pdata.porosity_initial_lava_flow,
        pdata.decay_depth_lava_flow,
        pdata.time_of_next_eruption_myr,
        pdata.eruption_interval_yr,
        pdata.iuse_eruption_interval,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
