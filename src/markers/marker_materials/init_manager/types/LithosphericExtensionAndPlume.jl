module LithosphericExtensionAndPlume

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import ..InitManager.InitStructs: LayerThickness
import ..InitManager.InitStructs: SevenLayerMaterialIDs
import ..InitManager.InitStructs: SeedCoordinates
import ..InitManager.InitStructs: PlumeGeometry
import ..InitManager.InitStructs: StrongZoneLimits
import ..InitManager.InitStructs: StrongZoneMaterialIDs
import ..SevenLayerEarthModel2D
import ..WeakSeed
import ..StrongRegions
import ..Plume: in_plume_region_check

function initialize!(
    model::ModelData;
    parameters::Union{Dict{String, Any}, Nothing}=nothing,
    material_domain_ids::Union{Dict{String, Int}, Nothing}=nothing
)::Nothing
    set_model_data!(model, parameters, material_domain_ids)
    initialize_material_ids!(model)
    return nothing
end

function initialize_material_ids!(model::ModelData)::Nothing
    earth_layering = model.geometry.parameters.earth_layering
    litho_strong_zones = model.geometry.parameters.litho_strong_zones
    plume = model.geometry.parameters.plume

    layer_thickness = LayerThickness(
        model.geometry.parameters.sticky_air_geometry.thick_air.value,
        earth_layering.thick_upper_crust.value,
        earth_layering.thick_lower_crust.value,
        earth_layering.thick_upper_lith.value,
        earth_layering.thick_middle_lith.value,
        earth_layering.thick_lower_lith.value
    )

    seed_coordinates = SeedCoordinates(
        model.geometry.parameters.weak_seed.x_seed.value,
        model.geometry.parameters.weak_seed.y_seed.value,
        model.geometry.parameters.weak_seed.w_seed.value
    )

    strong_zone_limits = StrongZoneLimits(
        litho_strong_zones.x_left_strong.value,
        litho_strong_zones.x_right_strong.value
    )

    plume_geometry = PlumeGeometry(
        plume.plume_radius.value,
        plume.plume_center_x.value,
        plume.plume_center_y.value,
        plume.plume_head_thick.value
    )

    matid_domains_dict = model.materials.dicts.matid_domains

    seven_layer_material_ids = SevenLayerMaterialIDs(
        matid_domains_dict["Atmosphere"],
        matid_domains_dict["UpperContinentalCrust"],
        matid_domains_dict["LowerContinentalCrust"],
        matid_domains_dict["UpperMantleLithosphere"],
        matid_domains_dict["MiddleMantleLithosphere"],
        matid_domains_dict["LowerMantleLithosphere"],
        matid_domains_dict["Asthenosphere"]
    )

    strong_zone_material_ids = StrongZoneMaterialIDs(
        matid_domains_dict["UpperContinentalCrustStrongZone"],
        matid_domains_dict["LowerContinentalCrustStrongZone"],
        matid_domains_dict["LithosphericMantleStrongZone"]
    )

    matid_weak_seed_mantle = matid_domains_dict["WeakSeedMantle"]
    matid_plume = matid_domains_dict["MantlePlume"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = define_material(
            x_marker, y_marker,
            layer_thickness, seven_layer_material_ids,
            seed_coordinates, matid_weak_seed_mantle,
            strong_zone_limits, strong_zone_material_ids,
            plume_geometry, matid_plume
        )

        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_material(
    x_marker::Float64,
    y_marker::Float64,
    layer_thickness::LayerThickness,
    seven_layer_material_ids::SevenLayerMaterialIDs,
    seed_coordinates::SeedCoordinates,
    matid_weak_seed_mantle::Int16,
    strong_zone_limits::StrongZoneLimits,
    strong_zone_material_ids::StrongZoneMaterialIDs,
    plume_geometry::PlumeGeometry,
    matid_plume::Int16;
    iuse_weak_seed::Int = 0,
    iuse_plume_head::Int = 1
)::Int16
    matid = SevenLayerEarthModel2D.define_material(
        y_marker,
        layer_thickness,
        seven_layer_material_ids
    )

    if iuse_weak_seed == 1
        matid = WeakSeed.define_material(
            matid,
            x_marker,
            y_marker,
            seed_coordinates,
            matid_weak_seed_mantle
        )
    else
        matid = StrongRegions.define_material(
            x_marker,
            y_marker,
            matid,
            strong_zone_limits,
            strong_zone_material_ids,
            layer_thickness
        )
    end

    if iuse_plume_head == 1
        if in_plume_region_check(x_marker, y_marker, plume_geometry)
            matid = matid_plume
        end
    end

    return matid
end

end # module LithosphericExtensionAndPlume 