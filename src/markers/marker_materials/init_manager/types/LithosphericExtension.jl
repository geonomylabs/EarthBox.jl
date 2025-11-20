module LithosphericExtension

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import ..InitManager.InitStructs: LayerThickness, SevenLayerMaterialIDs
import ..InitManager.InitStructs: StrongZoneLimits, StrongZoneMaterialIDs
import ..InitManager.InitStructs: SeedCoordinates, PlumeGeometry, WeakFaultGeometry
import ....Layers: calc_depths_for_layers
import ..SevenLayerEarthModel2D: define_material as define_seven_layer_material
import ..StrongRegions: define_material as define_strong_region_material
import ..WeakSeed: define_material as define_weak_seed_material
import ..Plume: in_plume_region_check
import ..LithosphericExtensionWeakFault: define_weak_fault_region

function initialize!(
    model::ModelData;
    parameters::Union{Dict{String, Any}, Nothing}=nothing,
    material_domain_ids::Union{Dict{String, Int}, Nothing}=nothing
)::Nothing
    set_model_data!(model, parameters, material_domain_ids)
    initialize_material_ids!(model)
    return nothing
end

function initialize_material_ids!(model::ModelData)
    (
        layer_thickness, strong_zone_limits, 
        seven_layer_material_ids, strong_zone_material_ids
    ) = get_lsz_material_model_inputs(model)
    matid_plume, plume_geometry = get_plume_model_inputs(model)
    matid_weak_seed_mantle, seed_coordinates = get_weak_seed_model_inputs(model)
    (
        matid_weak_crustal_fault, matid_weak_mantle_fault, 
        weak_fault_geometry
    ) = get_weak_fault_inputs(model)
    
    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    iuse_strong_zones = get_iuse_strong_zones(model)
    iuse_weak_fault = get_iuse_weak_fault(model)
    iuse_weak_seed = get_iuse_weak_seed(model)
    iuse_plume = get_iuse_plume(model)

    print_info("Material options:", level=1)
    print_info("iuse_strong_zones: $iuse_strong_zones", level=2)
    print_info("iuse_weak_fault: $iuse_weak_fault", level=2)
    print_info("iuse_weak_seed: $iuse_weak_seed", level=2)
    print_info("iuse_plume: $iuse_plume", level=2)

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        
        matid = define_material(
            x_marker, 
            y_marker,
            strong_zone_limits, 
            strong_zone_material_ids,
            layer_thickness, 
            seven_layer_material_ids,
            matid_plume, 
            plume_geometry,
            matid_weak_seed_mantle, 
            seed_coordinates,
            matid_weak_crustal_fault, 
            matid_weak_mantle_fault, 
            weak_fault_geometry,
            iuse_strong_zones = iuse_strong_zones,
            iuse_weak_fault = iuse_weak_fault,
            iuse_weak_seed = iuse_weak_seed,
            iuse_plume = iuse_plume
        )
        
        model.markers.arrays.material.marker_matid.array[imarker] = Int16(matid)
    end
end

function get_plume_model_inputs(model::ModelData)::Tuple{Int16, PlumeGeometry}
    matid_domains_dict = model.materials.dicts.matid_domains
    matid_plume = matid_domains_dict["MantlePlume"]

    plume = model.geometry.parameters.plume
    plume_geometry = PlumeGeometry(
        plume.plume_radius.value,
        plume.plume_center_x.value,
        plume.plume_center_y.value,
        plume.plume_head_thick.value
    )

    return matid_plume, plume_geometry
end

function get_weak_seed_model_inputs(model::ModelData)::Tuple{Int16, SeedCoordinates}
    matid_domains_dict = model.materials.dicts.matid_domains
    matid_weak_seed_mantle = matid_domains_dict["WeakSeedMantle"]

    weak_seed = model.geometry.parameters.weak_seed
    seed_coordinates = SeedCoordinates(
        weak_seed.x_seed.value,
        weak_seed.y_seed.value,
        weak_seed.w_seed.value
    )
    return matid_weak_seed_mantle, seed_coordinates
end

function get_weak_fault_inputs(model::ModelData)::Tuple{Int16, Int16, WeakFaultGeometry}
    weak_fault = model.geometry.parameters.weak_fault
    weak_fault_geometry = WeakFaultGeometry(
        weak_fault.fault_dip_degrees.value,
        weak_fault.fault_thickness.value,
        weak_fault.x_initial_fault.value,
        weak_fault.fault_height.value
    )
    
    matid_domains_dict = model.materials.dicts.matid_domains
    matid_weak_crustal_fault = matid_domains_dict["WeakCrustalFaultZone"]
    matid_weak_mantle_fault = matid_domains_dict["WeakMantleFaultZone"]
    
    return matid_weak_crustal_fault, matid_weak_mantle_fault, weak_fault_geometry
end

function get_lsz_material_model_inputs(
    model::ModelData
)::Tuple{LayerThickness, StrongZoneLimits, SevenLayerMaterialIDs, StrongZoneMaterialIDs}
    layer_thickness = get_layer_thickness(model)
    strong_zone_limits = get_strong_zone_limits(model)
    seven_layer_material_ids = get_seven_layer_material_ids(model)
    strong_zone_material_ids = get_strong_zone_material_ids(model)
    
    return (layer_thickness, strong_zone_limits, seven_layer_material_ids, strong_zone_material_ids)
end

function get_iuse_strong_zones(model::ModelData)::Int
    return model.geometry.parameters.litho_strong_zones.iuse_strong_zones.value
end

function get_iuse_weak_fault(model::ModelData)::Int
    return model.geometry.parameters.weak_fault.iuse_weak_fault.value
end

function get_iuse_weak_seed(model::ModelData)::Int
    return model.geometry.parameters.weak_seed.iuse_weak_seed.value
end

function get_iuse_plume(model::ModelData)::Int
    return model.geometry.parameters.plume.iuse_plume.value
end

function get_layer_thickness(model::ModelData)::LayerThickness
    earth_layering = model.geometry.parameters.earth_layering
    return LayerThickness(
        model.geometry.parameters.sticky_air_geometry.thick_air.value,
        earth_layering.thick_upper_crust.value,
        earth_layering.thick_lower_crust.value,
        earth_layering.thick_upper_lith.value,
        earth_layering.thick_middle_lith.value,
        earth_layering.thick_lower_lith.value
    )
end

function get_strong_zone_limits(model::ModelData)::StrongZoneLimits
    litho_strong_zones = model.geometry.parameters.litho_strong_zones
    return StrongZoneLimits(
        litho_strong_zones.x_left_strong.value,
        litho_strong_zones.x_right_strong.value
    )
end

function get_seven_layer_material_ids(model::ModelData)::SevenLayerMaterialIDs
    matid_domains_dict = model.materials.dicts.matid_domains
    return SevenLayerMaterialIDs(
        matid_domains_dict["Atmosphere"],
        matid_domains_dict["UpperContinentalCrust"],
        matid_domains_dict["LowerContinentalCrust"],
        matid_domains_dict["UpperMantleLithosphere"],
        matid_domains_dict["MiddleMantleLithosphere"],
        matid_domains_dict["LowerMantleLithosphere"],
        matid_domains_dict["Asthenosphere"]
    )
end

function get_strong_zone_material_ids(model::ModelData)::StrongZoneMaterialIDs
    matid_domains_dict = model.materials.dicts.matid_domains
    return StrongZoneMaterialIDs(
        matid_domains_dict["UpperContinentalCrustStrongZone"],
        matid_domains_dict["LowerContinentalCrustStrongZone"],
        matid_domains_dict["LithosphericMantleStrongZone"]
    )
end

function define_material(
    x_marker::Float64, 
    y_marker::Float64,
    strong_zone_limits::StrongZoneLimits,
    strong_zone_material_ids::StrongZoneMaterialIDs,
    layer_thickness::LayerThickness,
    seven_layer_material_ids::SevenLayerMaterialIDs,
    matid_plume::Int16,
    plume_geometry::PlumeGeometry,
    matid_weak_seed_mantle::Int16,
    seed_coordinates::SeedCoordinates,
    matid_weak_crustal_fault::Int16,
    matid_weak_mantle_fault::Int16,
    weak_fault_geometry::WeakFaultGeometry;
    iuse_strong_zones::Int = 1,
    iuse_weak_fault::Int = 0,
    iuse_weak_seed::Int = 0,
    iuse_plume::Int = 0
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

    # Apply plume if enabled
    if iuse_plume == 1
        if in_plume_region_check(x_marker, y_marker, plume_geometry)
            matid = matid_plume
        end
    end

    # Apply weak seed if enabled
    if iuse_weak_seed == 1
        matid = define_weak_seed_material(
            matid, x_marker, y_marker,
            seed_coordinates, matid_weak_seed_mantle
        )
    end

    # Apply weak fault if enabled
    if iuse_weak_fault == 1
        matid = define_weak_fault_region(
            matid, x_marker, y_marker,
            weak_fault_geometry, matid_weak_crustal_fault,
            matid_weak_mantle_fault, layer_thickness
        )
    end

    return matid
end

end # module LithosphericExtension 