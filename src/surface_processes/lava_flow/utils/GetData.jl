module GetData

import EarthBox.ModelDataContainer: ModelData
import EarthBox: MarkerSwarm
import EarthBox.Markers.MarkerMaterials: MaterialGroupIDs
import EarthBox.FindShallowest: find_shallowest_gabbroic_partially_molten_or_magma_marker

""" Get the extrusion parameters from the model data.

# Returns
- `y_sealevel`: The y-coordinate of the sea level (meters)
- `characteristic_flow_length_subaerial`: The characteristic flow length for subaerial flows (meters)
- `characteristic_flow_length_submarine`: The characteristic flow length for submarine flows (meters)
- `residual_lava_thickness_subaerial`: The residual lava thickness for subaerial flows (meters)
- `residual_laval_thickness_submarine`: The residual lava thickness for submarine flows (meters)
- `use_random_eruption_location`: A flag to use a random eruption location
- `use_normal_eruption_location`: A flag to use normal eruption location
- `decimation_factor`: The decimation factor for the extrusion model
"""
function get_extrusion_parameters(
    model::ModelData
)::Tuple{Float64, Float64, Float64, Float64, Float64, Bool, Bool, Int64}
    
    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    extrusion_params = model.melting.parameters.extrusion

    iuse_random_eruption_location = 
        extrusion_params.iuse_random_eruption_location.value
    use_random_eruption_location = iuse_random_eruption_location != 0

    iuse_normal_eruption_location = 
        extrusion_params.iuse_normal_eruption_location.value
    use_normal_eruption_location = iuse_normal_eruption_location != 0

    decimation_factor = extrusion_params.decimation_factor.value
    characteristic_flow_length_subaerial = 
        extrusion_params.characteristic_flow_length_subaerial.value
    characteristic_flow_length_submarine = 
        extrusion_params.characteristic_flow_length_submarine.value
    residual_lava_thickness_subaerial = 
        extrusion_params.residual_lava_thickness_subaerial.value
    residual_laval_thickness_submarine = 
        extrusion_params.residual_lava_thickness_submarine.value

    (
        characteristic_flow_length_subaerial,
        characteristic_flow_length_submarine,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine
    ) = modify_flow_parameters_for_magma_flush(
        model,
        characteristic_flow_length_subaerial,
        characteristic_flow_length_submarine,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine
    )

    return (
        y_sealevel,
        characteristic_flow_length_subaerial,
        characteristic_flow_length_submarine,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine,
        use_random_eruption_location,
        use_normal_eruption_location,
        decimation_factor
        )
end

""" Get the extrusion location parameters from the model data.

This function calculates the minimum x-coordinate of the eruption domain and the width of 
the eruption domain. The eruption domain is centered around the shallowest gabbroic 
partially molten or magma marker in the drainage basin.

# Arguments
- `model`: The model data structure
- `idrainage_basin`: The drainage basin index (default: -1)

# Returns
- `xmid_molten_domain`: The x-coordinate of the middle of the molten domain (meters)
- `width_eruption_domain`: The width of the eruption domain (meters)
- `eruption_location_x_min`: The minimum x-coordinate of the eruption location (meters)
"""
function get_extrusion_location_parameters(
    model::ModelData, 
    idrainage_basin::Int=-1
)::Tuple{Float64, Float64, Float64}
    
    xstart = model.melting.arrays.extraction.xstart_drainage.array[idrainage_basin]
    xend = model.melting.arrays.extraction.xend_drainage.array[idrainage_basin]
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    width_eruption_domain_fixed = calculate_width_eruption_domain_fixed(model)

    width_molten_domain = abs(
        model.melting.arrays.extraction.width_molten_zones.array[idrainage_basin])
    xmid_molten_zone = abs(
        model.melting.arrays.extraction.xmid_molten_zones.array[idrainage_basin])
    xshallow_partial_melt_avg = 
        model.melting.arrays.extraction.avg_shallow_partial_melt_xcoors.array[idrainage_basin]

    molten_gabbro_ids = MaterialGroupIDs.get_molten_gabbro_ids(model)

    (
        _imarker_magma_shallow, x_magma_shallow, _y_magma_shallow
    ) = find_shallowest_gabbroic_partially_molten_or_magma_marker(
            marker_x, marker_y, marker_matid,
            molten_gabbro_ids,
            xstart=xstart, xend=xend
        )

    (
        _xmid_gabbro_glacier, _ytop_gabbro_glacier, width_gabbro_glacier
    ) = calculate_dimensions_of_gabbro_glacier_domain(model, xstart=xstart, xend=xend)

    width_eruption_domain = if width_eruption_domain_fixed > 0.0
        width_eruption_domain_fixed
    elseif width_gabbro_glacier > 0.0 && abs(width_gabbro_glacier) < 1e32
        min(width_gabbro_glacier, width_molten_domain)
    else
        width_molten_domain
    end

    use_x_magma_shallow = false
    use_xmid_molten_zone = false
    eruption_location_x_min = if use_x_magma_shallow
        x_magma_shallow - width_eruption_domain/2.0
    elseif use_xmid_molten_zone
        xmid_molten_zone - width_eruption_domain/2.0
    else
        xshallow_partial_melt_avg - width_eruption_domain/2.0
    end

    (
        x_magma_shallow,
        width_eruption_domain,
        eruption_location_x_min
    ) = modify_flow_location_for_magma_flush(
        model,
        x_magma_shallow,
        width_eruption_domain,
        eruption_location_x_min
    )

    print_info = false
    if print_info
        print_extrusion_location_info(
            idrainage_basin, xstart, xend, _imarker_magma_shallow,
            x_magma_shallow, _y_magma_shallow, _xmid_gabbro_glacier,
            _ytop_gabbro_glacier, width_gabbro_glacier, width_molten_domain
        )
    end

    return (
        x_magma_shallow,
        width_eruption_domain,
        eruption_location_x_min
        )
end

""" Calculate width of eruption domain based on calculated characteristic magmatic crust height.

# Returns
- `width_eruption_domain_fixed`: The fixed width of the eruption domain
"""
function calculate_width_eruption_domain_fixed(model::ModelData)::Float64
    extrusion = model.melting.parameters.extrusion

    characteristic_magmatic_crust_height = 
        extrusion.characteristic_magmatic_crust_height.value

    width_eruption_domain_fixed_min = extrusion.width_eruption_domain_fixed.value
    width_eruption_domain_fixed_max = extrusion.width_eruption_domain_fixed_max.value
    if width_eruption_domain_fixed_max < width_eruption_domain_fixed_min
        width_eruption_domain_fixed_max = width_eruption_domain_fixed_min
    end

    characteristic_magmatic_crust_height_min = 
        extrusion.characteristic_magmatic_crust_height_min.value
    characteristic_magmatic_crust_height_max = 
        extrusion.characteristic_magmatic_crust_height_max.value

    width_eruption_domain_fixed = if characteristic_magmatic_crust_height <= 
        characteristic_magmatic_crust_height_min
        width_eruption_domain_fixed_min
    elseif characteristic_magmatic_crust_height >= characteristic_magmatic_crust_height_max
        width_eruption_domain_fixed_max
    else
        width_eruption_domain_fixed_min + 
            (width_eruption_domain_fixed_max - width_eruption_domain_fixed_min) /
            (characteristic_magmatic_crust_height_max - characteristic_magmatic_crust_height_min) *
            (characteristic_magmatic_crust_height - characteristic_magmatic_crust_height_min)
    end
    return width_eruption_domain_fixed
end

function print_extrusion_location_info(
    idrainage_basin::Int, 
    xstart::Float64, xend::Float64,
    _imarker_magma_shallow::Int, 
    x_magma_shallow::Float64, 
    _y_magma_shallow::Float64,
    _xmid_gabbro_glacier::Float64, 
    _ytop_gabbro_glacier::Float64,
    width_gabbro_glacier::Float64, 
    width_molten_domain::Float64
)::Nothing
    println(">> idrainage_basin: ", idrainage_basin)
    println(">> xstart: ", xstart)
    println(">> xend: ", xend)
    println(">> _imarker_magma_shallow: ", _imarker_magma_shallow)
    println(">> x_magma_shallow: ", x_magma_shallow)
    println(">> _y_magma_shallow: ", _y_magma_shallow)
    println(">> _xmid_gabbro_glacier: ", _xmid_gabbro_glacier)
    println(">> _ytop_gabbro_glacier: ", _ytop_gabbro_glacier)
    println(">> width_gabbro_glacier: ", width_gabbro_glacier)
    println(">> width_molten_domain: ", width_molten_domain)
    return nothing
end

""" Modify the flow parameters for magma flush.

# Returns
- Tuple of modified flow parameters
"""
function modify_flow_parameters_for_magma_flush(
    model::ModelData,
    characteristic_flow_length_subaerial::Float64,
    characteristic_flow_length_submarine::Float64,
    residual_lava_thickness_subaerial::Float64,
    residual_laval_thickness_submarine::Float64;
    flush_flow_length::Float64=400_000.0,
    flush_flow_thickness::Float64=30.0
)::Tuple{Float64, Float64, Float64, Float64}
    
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    initial_magma_flush_steps = 
        model.melting.parameters.extrusion.initial_magma_flush_steps.value
    
    if initial_magma_flush_steps > 0 && ntimestep <= initial_magma_flush_steps
        characteristic_flow_length_subaerial = flush_flow_length
        characteristic_flow_length_submarine = flush_flow_length
        residual_lava_thickness_subaerial = flush_flow_thickness
        residual_laval_thickness_submarine = flush_flow_thickness
    end
    
    return (
        characteristic_flow_length_subaerial,
        characteristic_flow_length_submarine,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine
        )
end

""" Modify the eruption location parameters for magma flush.

# Returns
- Tuple of modified location parameters
"""
function modify_flow_location_for_magma_flush(
    model::ModelData, 
    x_magma_shallow::Float64,
    width_eruption_domain::Float64, 
    eruption_location_x_min::Float64
)::Tuple{Float64,Float64, Float64}
    
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    initial_magma_flush_steps = 
        model.melting.parameters.extrusion.initial_magma_flush_steps.value
    
    if initial_magma_flush_steps > 0 && ntimestep <= initial_magma_flush_steps
        x_magma_shallow = -1e39
        width_eruption_domain = 0.0
        eruption_location_x_min = -1e39
    end
    
    return (
        x_magma_shallow,
        width_eruption_domain,
        eruption_location_x_min
        )
end

""" Check dimensions and print them.
"""
function check_dimensions(
    x_magma_shallow::Float64, y_magma_shallow::Float64,
    xmid_gabbro_glacier::Float64, ytop_gabbro_glacier::Float64,
    width_gabbro_glacier::Float64, width_molten_domain::Float64)::Nothing
    
    println("Check dimensions...")
    println(">> x_magma_shallow: ", x_magma_shallow)
    println(">> y_magma_shallow: ", y_magma_shallow)
    println(">> xmid_gabbro_glacier: ", xmid_gabbro_glacier)
    println(">> ytop_gabbro_glacier: ", ytop_gabbro_glacier)
    println(">> width_gabbro_glacier: ", width_gabbro_glacier)
    println(">> width_molten_domain: ", width_molten_domain)
    return nothing
end 

""" Calculate dimensions of gabbro glacier zone.

# Arguments
- model::ModelData
    - Model data structure.
- xstart::Float64
    - Starting x-coordinate of drainage basin.
- xend::Float64
    - Ending x-coordinate of drainage basin.

# Returns
- xmid_gabbro_glacier::Float64
    - X-location of the middle of molten zone in meters.
- ytop_gabbro_glacier::Float64
    - Minimum y-location of the molten zone in meters.
- width_gabbro_glacier::Float64
    - Width of molten zone in meters.
"""
function calculate_dimensions_of_gabbro_glacier_domain(
    model::ModelData;
    xstart::Float64 = -1e39,
    xend::Float64 = 1e39
)::Tuple{Float64, Float64, Float64}
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    matid_types = model.materials.dicts.matid_types

    matid_solidified_layered_gabbro_partially_molten = matid_types[
        "SolidifiedLayeredGabbroPartiallyMolten"][1]
    target_matids = [matid_solidified_layered_gabbro_partially_molten]

    (
        xmid_gabbro_glacier, ytop_gabbro_glacier, width_gabbro_glacier
    ) = MarkerSwarm.calculate_dimensions_general(
        marker_x, marker_y, marker_matid,
        target_matids,
        xstart, xend
    )
    return xmid_gabbro_glacier, ytop_gabbro_glacier, width_gabbro_glacier
end


end # module