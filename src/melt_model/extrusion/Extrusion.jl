module Extrusion

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_extrusion::Symbol
    iuse_random_eruption_location::Symbol
    iuse_normal_eruption_location::Symbol
    iuse_eruption_interval::Symbol
    extrusion_volume_factor::Symbol
    extrusion_volume_factor_max::Symbol
    characteristic_magmatic_crust_height_min::Symbol
    characteristic_magmatic_crust_height_max::Symbol
    width_eruption_domain_fixed::Symbol
    width_eruption_domain_fixed_max::Symbol
    characteristic_flow_length_subaerial::Symbol
    characteristic_flow_length_submarine::Symbol
    residual_lava_thickness_subaerial::Symbol
    residual_lava_thickness_submarine::Symbol
    decimation_factor::Symbol
    initial_magma_flush_steps::Symbol
    magma_flush_factor::Symbol
    porosity_initial_lava_flow::Symbol
    decay_depth_lava_flow::Symbol
    eruption_interval_yr::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize melt extrusion model parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing the 
    model parameters and arrays.

# Keyword Arguments
- `$(PDATA.iuse_extrusion.name)::Int64`
    - $(PDATA.iuse_extrusion.description)
- `$(PDATA.iuse_random_eruption_location.name)::Int64`
    - $(PDATA.iuse_random_eruption_location.description)
- `$(PDATA.iuse_normal_eruption_location.name)::Int64`
    - $(PDATA.iuse_normal_eruption_location.description)
- `$(PDATA.iuse_eruption_interval.name)::Int64`
    - $(PDATA.iuse_eruption_interval.description)
- `$(PDATA.extrusion_volume_factor.name)::Float64`
    - $(PDATA.extrusion_volume_factor.description)
- `$(PDATA.extrusion_volume_factor_max.name)::Float64`
    - $(PDATA.extrusion_volume_factor_max.description)
- `$(PDATA.characteristic_magmatic_crust_height_min.name)::Float64`
    - $(PDATA.characteristic_magmatic_crust_height_min.description)
- `$(PDATA.characteristic_magmatic_crust_height_max.name)::Float64`
    - $(PDATA.characteristic_magmatic_crust_height_max.description)
- `$(PDATA.width_eruption_domain_fixed.name)::Float64`
    - $(PDATA.width_eruption_domain_fixed.description)
- `$(PDATA.width_eruption_domain_fixed_max.name)::Float64`
    - $(PDATA.width_eruption_domain_fixed_max.description)
- `$(PDATA.residual_lava_thickness_subaerial.name)::Float64`
    - $(PDATA.residual_lava_thickness_subaerial.description)
- `$(PDATA.characteristic_flow_length_submarine.name)::Float64`
    - $(PDATA.characteristic_flow_length_submarine.description)
- `$(PDATA.characteristic_flow_length_subaerial.name)::Float64`
    - $(PDATA.characteristic_flow_length_subaerial.description)
- `$(PDATA.residual_lava_thickness_submarine.name)::Float64`
    - $(PDATA.residual_lava_thickness_submarine.description)
- `$(PDATA.decimation_factor.name)::Int64`
    - $(PDATA.decimation_factor.description)
- `$(PDATA.initial_magma_flush_steps.name)::Int64`
    - $(PDATA.initial_magma_flush_steps.description)
- `$(PDATA.magma_flush_factor.name)::Float64`
    - $(PDATA.magma_flush_factor.description)
- `$(PDATA.porosity_initial_lava_flow.name)::Float64`
    - $(PDATA.porosity_initial_lava_flow.description)
- `$(PDATA.decay_depth_lava_flow.name)::Float64`
    - $(PDATA.decay_depth_lava_flow.description)
- `$(PDATA.eruption_interval_yr.name)::Float64`
    - $(PDATA.eruption_interval_yr.description)

"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

""" Set extrusion thickness at topography grid node to zero.

Extrusion thickness is set to zero at the beginning of each time step.

# Arguments
- `gridt::Array{Float64,2}`: Topography grid array. This function updates 
  gridt[7,:], which stores extrusion thickness (m) at topography node.
"""
function zero_out_extrusion_thickness_for_timestep(gridt::Array{Float64,2})
    toponum = size(gridt, 2)
    for j in 1:toponum
        gridt[7, j] = 0.0
    end
end

""" Update extrusion thickness at topography grid node.

# Arguments
- `total_lava_thickness::Array{Float64,1}`: Total lava thickness (m) at 
  topography grid node that was emplaced during time step.
- `gridt::Array{Float64,2}`: Topography grid array. This function updates 
  gridt[7,:], which stores extrusion thickness (m) at topography node.
"""
function update_extrusion_thickness(
    total_lava_thickness::Array{Float64,1},
    gridt::Array{Float64,2}
)::Nothing
    toponum = size(gridt, 2)
    for j in 1:toponum
        gridt[7, j] = gridt[7, j] + total_lava_thickness[j]
    end
    return nothing
end

""" Calculate thickness of extruded material at topography grid node.

This is an old approach that is no longer used.

# Arguments
- `model::ModelData`: Model data container containing topography and melting 
  parameters.
"""
function calculate_extrusion_thickness(model::ModelData)
    dx_topo = model.topography.parameters.topo_grid.dx_topo.value
    toponum = model.topography.parameters.topo_grid.toponum.value
    gridt = model.topography.arrays.gridt.array

    extrusion_volume = model.melting.parameters.extraction.ext_vol.value
    xmid_mol = model.melting.parameters.extraction.xmid_mol.value
    width_mol = model.melting.parameters.extraction.width_mol.value

    width_ext = width_mol/2.0
    xmin_ext = xmid_mol - width_ext/2.0
    xmax_ext = xmid_mol + width_ext/2.0
    icount = 0

    # Count nodes in extrusion window
    for j in 1:toponum
        xtopo = gridt[1, j]
        if xmin_ext < xtopo < xmax_ext
            icount += 1
        end
    end

    # Calculate thickness at each node
    for j in 1:toponum
        xtopo = gridt[1, j]
        if xmin_ext < xtopo < xmax_ext
            gridt[7, j] = extrusion_volume/icount/dx_topo
        end
    end
end

end # module 