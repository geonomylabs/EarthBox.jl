module ContractionViscoelastic

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Markers.MarkerMaterials.InitManager.InitStructs:
    LayerThickness, SevenLayerMaterialIDs, StrongZoneLimits, StrongZoneMaterialIDs
import EarthBox.Markers.MarkerMaterials.InitManager.LithosphericExtension:
    get_lsz_material_model_inputs
import EarthBox.Markers.MarkerMaterials.InitManager.SevenLayerEarthModel2D: define_material as define_seven_layer_material
import EarthBox.Markers.MarkerMaterials.InitManager.StrongRegions: define_material as define_strong_region_material
import ..Coordinates: MarkerRecycleArrays
import ..Coordinates: reset_recycled_marker_coordinates
import ..ResetMarkers: SubsurfaceProps
import ..ResetMarkers: make_recycle_property_arrays
import ..ResetMarkers: reset_recycled_marker_properties!
import ..PressureReset: calculate_reset_pressure
import ..ContractionRecycleLocations: calculate_recycled_locations
import ..InitializeMarkerRecycling: initialize_recycling
import ..CheckDicts: check_domain_dict

function recycle_markers_symmetric!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    check_required_types_and_domains_fracture_zone(model)
    execute_recycle_steps(model, inside_flags; use_asymmetric_contraction=false)
    return nothing
end

function recycle_markers_asymmetric!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    check_required_types_and_domains_fracture_zone(model)
    execute_recycle_steps(model, inside_flags; use_asymmetric_contraction=true)
    return nothing
end

function check_required_types_and_domains_fracture_zone(model::ModelData)::Nothing
    matid_types = model.materials.dicts.matid_types
    type_dict = Dict(
        "StickyAir" => matid_types["StickyAir"][1],
        "StickyWater" => matid_types["StickyWater"][1]
    )
    check_domain_dict(type_dict, "type", "contractional")

    domains = model.materials.dicts.matid_domains
    domain_dict = Dict(
        "Asthenosphere" => domains["Asthenosphere"],
        "UpperContinentalCrust" => domains["UpperContinentalCrust"],
        "LowerContinentalCrust" => domains["LowerContinentalCrust"],
        "UpperContinentalCrustStrongZone" => domains["UpperContinentalCrustStrongZone"],
        "LowerContinentalCrustStrongZone" => domains["LowerContinentalCrustStrongZone"],
        "UpperMantleLithosphere" => domains["UpperMantleLithosphere"],
        "MiddleMantleLithosphere" => domains["MiddleMantleLithosphere"],
        "LowerMantleLithosphere" => domains["LowerMantleLithosphere"],
        "LithosphericMantleStrongZone" => domains["LithosphericMantleStrongZone"]
    )
    check_domain_dict(domain_dict, "domain", "contractional")
    return nothing
end

function execute_recycle_steps(
    model::ModelData,
    inside_flags::Vector{Int8};
    use_asymmetric_contraction::Bool = false
)::Nothing
    (
        outside_indices, nrecycle_air, nrecycle_subsurface
    ) = initialize_recycling(model, inside_flags)
    recycled_coordinates = calculate_recycled_locations(model, nrecycle_air, nrecycle_subsurface, use_asymmetric_contraction)
    marker_recycling = reset_recycled_marker_coordinates(model, recycled_coordinates, outside_indices)
    subsurface_props = get_subsurface_reset_properties_viscoelastic_contraction(model, marker_recycling)
    reset_recycled_marker_properties!(model, marker_recycling, subsurface_props)
    return nothing
end

function get_subsurface_reset_properties_viscoelastic_contraction(
    model::ModelData,
    marker_recycling::MarkerRecycleArrays
)::SubsurfaceProps
    recycle_indices = marker_recycling.recycle_indices
    recycle_flags = marker_recycling.recycle_flags
    nrecycle_total = length(recycle_indices)
    (
        matids_recycle, pressure_recycle, temperature_recycle
    ) = make_recycle_property_arrays(nrecycle_total)

    marker_arrays = model.markers.arrays
    location = marker_arrays.location
    marker_x = location.marker_x.array
    marker_y = location.marker_y.array

    (
        layer_thickness, strong_zone_limits,
        seven_layer_material_ids, strong_zone_material_ids
    ) = get_lsz_material_model_inputs(model)

    for irecycle in 1:nrecycle_total
        recycle_flag = recycle_flags[irecycle]
        if recycle_flag == 1
            imarker = recycle_indices[irecycle]
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            matids_recycle[irecycle] = define_material(
                x_marker, y_marker,
                strong_zone_limits, strong_zone_material_ids,
                layer_thickness, seven_layer_material_ids
            )
            pressure_recycle[irecycle] = calculate_reset_pressure(y_marker)
            temperature_recycle[irecycle] = 0.0
        end
    end

    return SubsurfaceProps(matids_recycle, pressure_recycle, temperature_recycle)
end

function define_material(
    x_marker::Float64, 
    y_marker::Float64,
    strong_zone_limits::StrongZoneLimits,
    strong_zone_material_ids::StrongZoneMaterialIDs,
    layer_thickness::LayerThickness,
    seven_layer_material_ids::SevenLayerMaterialIDs;
    iuse_strong_zones::Int = 1
)::Int16
    # Start with the seven layer earth model
    matid = define_seven_layer_material(y_marker, layer_thickness, seven_layer_material_ids)

    # Apply strong zones if enabled
    if iuse_strong_zones == 1
        matid = define_strong_region_material(
            x_marker, y_marker, matid,
            strong_zone_limits, strong_zone_material_ids,
            layer_thickness
        )
    end

    return matid
end

end # module 