module UpdateHydrothermalCirculation

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox: MathTools
import EarthBox.Markers.MarkerMaterials: MaterialGroupIDs
import EarthBox: SedimentWaterInterface
import EarthBox.FindShallowest: find_shallowest_gabbroic_partially_molten_or_magma_marker
import EarthBox.SedimentThickness: calculate_sediment_thickness_from_markers_no_basalt
import EarthBox.ConversionFuncs: kelvin_to_celsius

""" Calculate rock properties for all markers.

Rock properties are a function of temperature, pressure, melt fraction,
submud depth, sediment thickness, water depth and composition.

# Updated Arrays
- `model.markers.arrays.thermal.marker_kt.array`: Array of marker thermal conductivity in W/m/K
"""
function update_for_hydrothermal_opt!(model::ModelData)
    # Option flags
    iuse_plastic_strain_rate_for_hydrothermal_cooling = model.materials.parameters.hydrothermal.iuse_plastic_strain_rate_for_hydrothermal.value
    iuse_plastic_strain_for_hydrothermal_cooling      = model.materials.parameters.hydrothermal.iuse_plastic_strain_for_hydrothermal.value
    iuse_hydrothermal                                 = model.materials.parameters.hydrothermal.iuse_hydrothermal.value
    iuse_topo                                         = model.topography.parameters.model_options.iuse_topo.value
    iuse_melt_lens                                    = model.materials.parameters.hydrothermal.iuse_melt_lens.value
    # Hydrothermal parameters
    hydrothermal_decay_length                  = model.materials.parameters.hydrothermal.hydrothermal_decay_length.value
    hydrothermal_buffer_distance               = model.materials.parameters.hydrothermal.hydrothermal_buffer_distance.value
    hydrothermal_plastic_strain_rate_reference = model.materials.parameters.hydrothermal.hydrothermal_plastic_strain_rate_reference.value
    hydrothermal_plastic_strain_reference      = model.materials.parameters.hydrothermal.hydrothermal_plastic_strain_reference.value
    smoothing_factor                           = model.materials.parameters.hydrothermal.hydrothermal_smoothing_factor.value
    nusselt_number                             = model.materials.parameters.hydrothermal.hydrothermal_nusselt_number.value
    max_temperature                            = model.materials.parameters.hydrothermal.hydrothermal_max_temperature.value
    max_depth                                  = model.materials.parameters.hydrothermal.hydrothermal_max_depth.value
    sediment_thickness_threshold               = model.materials.parameters.hydrothermal.sediment_thickness_threshold.value
    y_sealevel                                 = model.topography.parameters.sealevel.y_sealevel.value
    # Marker and drainage basin parameters
    marknum          = model.markers.parameters.distribution.marknum.value
    ndrainage_basin  = model.melting.parameters.extraction.ndrainage_basin.value
    # Arrays
    gridt                      = model.topography.arrays.gridt.array
    marker_strain_rate_plastic = model.markers.arrays.strain.marker_strain_rate_plastic.array
    marker_plastic_strain      = model.markers.arrays.strain.marker_strain_plastic.array
    marker_x                   = model.markers.arrays.location.marker_x.array
    marker_y                   = model.markers.arrays.location.marker_y.array
    marker_kt                  = model.markers.arrays.thermal.marker_kt.array
    marker_matid               = model.markers.arrays.material.marker_matid.array
    marker_TK                  = model.markers.arrays.thermal.marker_TK.array
    xstart_drainage            = model.melting.arrays.extraction.xstart_drainage.array
    xend_drainage              = model.melting.arrays.extraction.xend_drainage.array
    # Material group IDs
    matid_types = model.materials.dicts.matid_types
    ids_gabbro     = MaterialGroupIDs.get_gabbro_ids_array(model)
    ids_magma      = MaterialGroupIDs.get_magma_ids_array(model)
    ids_basalt     = matid_types["SolidifiedBasalt"]
    ids_mantle     = MaterialGroupIDs.get_mantle_ids_array(model)
    ids_cont_crust = MaterialGroupIDs.get_continental_crust_ids_array(model)

    topo_grid_sediment_thickness = calculate_sediment_thickness_from_markers_no_basalt(model)
    topo_gridx = gridt[1,:]
    topo_gridy = gridt[2,:]
    for idrainage_basin in 1:ndrainage_basin
        xstart = xstart_drainage[idrainage_basin]
        xend = xend_drainage[idrainage_basin]
        (
            xmin_lens, xmax_lens, _x_melt_lens, 
            y_melt_lens, y_top_gabbro_glacier,
            conductivity_reduction_factor, imarker_shallow
        ) = get_melt_lens_parameters(model, idrainage_basin)

        if imarker_shallow == -999
            iuse_melt_lens = 0
        end
        xmin_molten_domain, xmax_molten_domain = get_x_limits_of_molten_zone(model, idrainage_basin)
        print_info = false
        if print_info
            print_molten_domain_info(
                idrainage_basin, xmin_lens, xmax_lens,
                xmin_molten_domain, xmax_molten_domain
            )
        end
        if use_hydrothermal_loop(iuse_hydrothermal, iuse_topo, xmin_molten_domain, xmax_molten_domain)
            Threads.@threads for imarker in 1:marknum
                x_marker = marker_x[imarker]
                if xstart < x_marker < xend
                    @inbounds begin
                        y_marker            = marker_y[imarker]
                        material_id         = marker_matid[imarker]
                        strain_plastic      = marker_plastic_strain[imarker]
                        strain_rate_plastic = marker_strain_rate_plastic[imarker]
                        temperature_kelvins = marker_TK[imarker]
                        conductivity        = marker_kt[imarker]
                    end
                    sediment_thickness = MathTools.linear_interp_bisection(
                        topo_gridx, topo_grid_sediment_thickness, x_marker)
                    y_mudline = MathTools.linear_interp_bisection(
                        topo_gridx, topo_gridy, x_marker)
                    if use_melt_lens_model(
                        iuse_melt_lens, x_marker, y_marker, y_mudline,
                        y_top_gabbro_glacier, xmin_lens, xmax_lens,
                        y_melt_lens
                    )
                        conductivity = conductivity/conductivity_reduction_factor
                        marker_kt[imarker] = conductivity
                    elseif use_crustal_hydrothermal_model(
                        ids_gabbro, ids_magma, ids_basalt, ids_mantle,
                        ids_cont_crust, material_id, y_sealevel, y_mudline,
                        y_marker, sediment_thickness, sediment_thickness_threshold
                    )
                        nusselt_number_calculated = 
                            limit_max_nusselt_number_using_distance_from_molten_zone(
                                x_marker, nusselt_number, xmin_molten_domain,
                                xmax_molten_domain, hydrothermal_decay_length,
                                hydrothermal_buffer_distance
                            )

                        if iuse_plastic_strain_rate_for_hydrothermal_cooling == 1
                            nusselt_number_calculated = 
                                limit_max_nusselt_number_using_plastic_strain(
                                    nusselt_number_calculated, strain_plastic,
                                    strain_rate_plastic,
                                    hydrothermal_plastic_strain_rate_reference,
                                    hydrothermal_plastic_strain_reference,
                                    iuse_plastic_strain_for_hydrothermal_cooling
                                )
                        end
                        conductivity = calculate_rock_props(
                            temperature_kelvins, conductivity,
                            y_marker, y_mudline, smoothing_factor,
                            nusselt_number_calculated, max_temperature,
                            max_depth
                        )

                        marker_kt[imarker] = conductivity
                    end
                end
            end
        end
    end
end

""" Calculate thermal conductivity for hydrothermal circulation model.

Conductivity is a function of Nusselt number, temperature and sub-mud depth.

conductivity = (
    conductivity 
    + conductivity
    *(nusselt_number - 1.0)
    *exp(smoothing_factor*(1 - temperature/max_temperature))
    *exp(smoothing_factor*(1 - y_submud/max_depth))
)

or

conductivity = (
    conductivity
    + conductivity
    *(nusselt_number - 1.0)
    *exp(smoothing_factor*(2 - temperature/max_temperature - y_submud/max_depth))
)

# Arguments
- `temperature_kelvins::Float64`: Marker temperature in Kelvin
- `conductivity::Float64`: Thermal conductivity of the rock without hydrothermal circulation (W/m/K)
- `y_location::Float64`: Y location of marker in meters
- `y_mudline::Float64`: Y location of mudline in meters
- `smoothing_factor::Float64`: Smoothing factor for conductivity used in hydrothermal circulation model
- `nusselt_number::Float64`: Nusselt number used in hydrothermal circulation model
- `max_temperature::Float64`: Maximum temperature (Celsius) using in hydrothermal circulation model
- `max_depth::Float64`: Maximum depth in meters used in hydrothermal circulation model

# Returns
- `conductivity::Float64`: Updated thermal conductivity of the rock (W/m/K) used to approximate the effects of hydrothermal circulation
"""
@inline function calculate_rock_props(
    temperature_kelvins::Float64,
    conductivity::Float64,
    y_location::Float64,
    y_mudline::Float64,
    smoothing_factor::Float64,
    nusselt_number::Float64,
    max_temperature::Float64,
    max_depth::Float64
)::Float64
    y_marker = y_location
    temperature_celsius = kelvin_to_celsius(temperature_kelvins)
    if temperature_celsius < max_temperature
        y_submud = y_marker - y_mudline
        nusselt_number = max(1.0, nusselt_number)
        if 0.0 < y_submud < max_depth
            term1 = conductivity * (nusselt_number - 1.0)
            inner_term = 2.0 - temperature_celsius / max_temperature - y_submud / max_depth
            term2 = exp(smoothing_factor * inner_term)
            conductivity += term1 * term2
        end
    end
    return conductivity
end 

""" Get x limits of molten zone.

# Returns
- `Tuple{Float64,Float64}`: Minimum and maximum x-coordinates of molten domain
"""
function get_x_limits_of_molten_zone(
    model::ModelData,
    idrainage_basin::Int
)::Tuple{Float64,Float64}
    width_molten_zones = 
        model.melting.arrays.extraction.width_molten_zones.array
    width_molten_domain = width_molten_zones[idrainage_basin]

    xmid_molten_zones = 
        model.melting.arrays.extraction.xmid_molten_zones.array
    xmid_molten_domain = xmid_molten_zones[idrainage_basin]

    if width_molten_domain > 0.0
        xmin_molten_domain = xmid_molten_domain - width_molten_domain/2
        xmax_molten_domain = xmid_molten_domain + width_molten_domain/2
    else
        xmin_molten_domain = -1e39
        xmax_molten_domain = -1e39
    end

    return xmin_molten_domain, xmax_molten_domain
end

""" Get melt lens parameters.

# Returns
- `Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64}`: Parameters for melt lens
"""
function get_melt_lens_parameters(
    model::ModelData,
    idrainage_basin::Int
)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64}
    gridt = model.topography.arrays.gridt.array

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    xstart = model.melting.arrays.extraction.xstart_drainage.array[idrainage_basin]
    xend = model.melting.arrays.extraction.xend_drainage.array[idrainage_basin]

    molten_layered_gabbro_ids = 
        MaterialGroupIDs.get_molten_layered_gabbro_ids(model)
    imarker_shallow, x_top_gabbro_glacier, y_top_gabbro_glacier = 
        find_shallowest_gabbroic_partially_molten_or_magma_marker(
            marker_x, marker_y, marker_matid,
            molten_layered_gabbro_ids,
            xstart=xstart, xend=xend
        )
    y_mudline = SedimentWaterInterface.get_depth(x_top_gabbro_glacier, gridt)
    y_submud_top_gabbro_glacier = y_top_gabbro_glacier - y_mudline

    x_melt_lens = x_top_gabbro_glacier
    y_melt_lens = y_mudline + (y_submud_top_gabbro_glacier)*0.25

    surface_to_melt_lens = y_melt_lens - y_mudline

    width_melt_lens = 2_500.0
    xmin_lens = x_top_gabbro_glacier - width_melt_lens/2
    xmax_lens = x_top_gabbro_glacier + width_melt_lens/2

    conductivity_reduction_factor = 10.0
    print_info = false
    if print_info
        print_melt_lens_info(
            x_top_gabbro_glacier, y_top_gabbro_glacier,
            xmin_lens, xmax_lens, y_melt_lens,
            surface_to_melt_lens, conductivity_reduction_factor
        )
    end
    return (
        xmin_lens, xmax_lens, x_melt_lens,
        y_melt_lens, y_top_gabbro_glacier,
        conductivity_reduction_factor,
        imarker_shallow
    )
end

""" Check if hydrothermal circulation should be used to modify properties.

# Returns
- `Bool`: Flag indicating whether hydrothermal circulation loop should be used
"""
function use_hydrothermal_loop(
    iuse_hydrothermal::Int,
    iuse_topo::Int,
    xmin_molten_domain::Float64,
    xmax_molten_domain::Float64
)::Bool
    check = false
    if iuse_topo == 1 && iuse_hydrothermal == 1
        check = true
    end
    if xmin_molten_domain < 0.0 && xmax_molten_domain < 0.0
        check = false
    end
    return check
end

""" Check if melt lens model should be used to modify properties.

The melt lens parameterization is used to increase thermal conductivity above and 
below the melt lens in order to approximate the effects of melt advection to a 
shallow melt lens followed by downward flow of fractionated crystal mush.

# Returns
- `Bool`: Flag indicating whether melt lens model should be used
"""
function use_melt_lens_model(
    iuse_melt_lens::Int,
    x_marker::Float64,
    y_marker::Float64,
    y_mudline::Float64,
    y_top_gabbro_glacier::Float64,
    xmin_lens::Float64,
    xmax_lens::Float64,
    y_melt_lens::Float64
)::Bool
    check = false
    if iuse_melt_lens == 1
        if 0.0 < y_melt_lens < 50_000.0 && xmin_lens >= 0.0
            if y_mudline < y_marker < y_top_gabbro_glacier
                if xmin_lens < x_marker < xmax_lens
                    check = true
                end
            end
        end
    end
    return check
end

""" Check if hydrothermal circulation should be used to modify properties.

# Returns
- `Bool`: Flag indicating whether hydrothermal circulation should be used
"""
function use_crustal_hydrothermal_model(
    ids_gabbro::Vector{Int16},
    ids_magma::Vector{Int16},
    ids_basalt::Vector{Int16},
    ids_mantle::Vector{Int16},
    ids_cont_crust::Vector{Int16},
    material_id::Int16,
    y_sealevel::Float64,
    y_mudline::Float64,
    y_marker::Float64,
    sediment_thickness::Float64,
    sediment_thickness_threshold::Float64
)::Bool
    check = false
    if material_id in ids_gabbro ||
       material_id in ids_magma ||
       material_id in ids_basalt ||
       material_id in ids_mantle ||
       material_id in ids_cont_crust
        check = true
        if y_sealevel > y_mudline ||
           y_marker < y_sealevel ||
           sediment_thickness > sediment_thickness_threshold
            check = false
        end
    end
    return check
end

""" Calculate hydrothermal Nusselt number based on distance from molten zone.

# Returns
- `Float64`: Modified Nusselt number
"""
@inline function limit_max_nusselt_number_using_distance_from_molten_zone(
    x_location::Float64,
    nusselt_number::Float64,
    xmin_molten_domain::Float64,
    xmax_molten_domain::Float64,
    hydrothermal_decay_length::Float64,
    hydrothermal_buffer_distance::Float64
)::Float64
    hydrothermal_distance_factor = calculate_hydrothermal_distance_factor(
        x_location, xmin_molten_domain, xmax_molten_domain,
        hydrothermal_decay_length, hydrothermal_buffer_distance
    )
    nusselt_number = max(1.0, nusselt_number*hydrothermal_distance_factor)
    return nusselt_number
end

""" Calculate distance from molten zone factor.

# Returns
- `Float64`: Factor to modify properties based on distance from molten zone
"""
@inline function calculate_hydrothermal_distance_factor(
    x_location::Float64,
    xmin_molten_domain::Float64,
    xmax_molten_domain::Float64,
    hydrothermal_decay_length::Float64,
    hydrothermal_buffer_distance::Float64
)::Float64
    xmin_molten_domain = xmin_molten_domain - hydrothermal_buffer_distance
    xmax_molten_domain = xmax_molten_domain + hydrothermal_buffer_distance
    distance_from_molten_zone = calculate_distance_from_molten_zone(
        x_location, xmin_molten_domain, xmax_molten_domain)
    factor = exp(-distance_from_molten_zone/hydrothermal_decay_length)
    return factor
end

""" Calculate distance from molten zone.

# Returns
- `Float64`: Distance from molten zone
"""
@inline function calculate_distance_from_molten_zone(
    x_location::Float64,
    xmin_molten_domain::Float64,
    xmax_molten_domain::Float64
)::Float64
    distance = 0.0
    if xmin_molten_domain < x_location < xmax_molten_domain
        distance = 0.0
    elseif x_location < xmin_molten_domain
        distance = xmin_molten_domain - x_location
    elseif x_location > xmax_molten_domain
        distance = x_location - xmax_molten_domain
    end
    return distance
end

""" Calculate Nusselt number for hydrothermal circulation using plastic strain rate.

# Returns
- `Float64`: Updated Nusselt number
"""
@inline function limit_max_nusselt_number_using_plastic_strain(
    nusselt_number_max::Float64,
    strain_plastic::Float64,
    strain_rate_plastic::Float64,
    hydrothermal_plastic_strain_rate_reference::Float64,
    strain_plastic_reference::Float64,
    iuse_plastic_strain::Int
)::Float64
    strain_rate_ratio = 
        strain_rate_plastic/hydrothermal_plastic_strain_rate_reference

    nusselt_number = calculate_nusselt_number_using_strain_ratio(
        strain_rate_ratio, nusselt_number_max)

    if iuse_plastic_strain == 1
        strain_ratio = strain_plastic/strain_plastic_reference
        nusselt_number_strain = calculate_nusselt_number_using_strain_ratio(
            strain_ratio, nusselt_number_max)
        nusselt_number = max(nusselt_number, nusselt_number_strain)
    end

    nusselt_number = max(1.0, nusselt_number)
    return nusselt_number
end

""" Calculate Nusselt number for hydrothermal circulation using a ratio.

# Returns
- `Float64`: Calculated Nusselt number
"""
@inline function calculate_nusselt_number_using_strain_ratio(
    strain_ratio::Float64,
    nusselt_number_max::Float64
)::Float64
    nusselt_number = 
        1.0 + (nusselt_number_max - 1.0)*(1.0 - exp(-strain_ratio))
    return nusselt_number
end

function print_molten_domain_info(
    idrainage_basin::Int,
    xmin_lens::Float64,
    xmax_lens::Float64,
    xmin_molten_domain::Float64,
    xmax_molten_domain::Float64
)
    print_info(">> Working on hydrothermal circulation for drainage basin $(idrainage_basin)", level=1)
    print_info(">> Molten x-limits for hydrothermal circulation:", level=2)
    print_info("xmin_melt_lens: $(xmin_lens)", level=2)
    print_info("xmax_melt_lens: $(xmax_lens)", level=2)
    print_info("xmin_molten_domain: $(xmin_molten_domain)", level=2)
    print_info("xmax_molten_domain: $(xmax_molten_domain)", level=2)
    print_info("", level=2)
end

function print_melt_lens_info(
    x_top_gabbro_glacier::Float64,
    y_top_gabbro_glacier::Float64,
    xmin_lens::Float64,
    xmax_lens::Float64,
    y_melt_lens::Float64,
    surface_to_melt_lens::Float64,
    conductivity_reduction_factor::Float64
)
    print_info("", level=1)
    print_info(">> Melt lens parameters:", level=2)
    print_info("x_top_gabbro_glacier: $(x_top_gabbro_glacier)", level=2)
    print_info("y_top_gabbro_glacier: $(y_top_gabbro_glacier)", level=2)
    print_info("xmin_lens: $(xmin_lens)", level=2)
    print_info("xmax_lens: $(xmax_lens)", level=2)
    print_info("y_melt_lens: $(y_melt_lens)", level=2)
    print_info("surface_to_melt_lens: $(surface_to_melt_lens)", level=2)
    print_info("conductivity_reduction_factor: $(conductivity_reduction_factor)", level=2)
    print_info("", level=1)
end

end # module 