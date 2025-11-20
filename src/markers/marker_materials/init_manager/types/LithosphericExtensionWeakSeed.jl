module LithosphericExtensionWeakSeed

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import ..InitManager.InitStructs: LayerThickness
import ..InitManager.InitStructs: SevenLayerMaterialIDs
import ..InitManager.InitStructs: SeedCoordinates
import ..SevenLayerEarthModel2D
import ..WeakSeed

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

    matid_weak_seed_mantle = matid_domains_dict["WeakSeedMantle"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = define_material(
            x_marker, y_marker, layer_thickness,
            seven_layer_material_ids, seed_coordinates,
            matid_weak_seed_mantle
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
    matid_weak_seed_mantle::Int16
)::Int16
    matid = SevenLayerEarthModel2D.define_material(
        y_marker,
        layer_thickness,
        seven_layer_material_ids
    )

    matid = WeakSeed.define_material(
        matid,
        x_marker,
        y_marker,
        seed_coordinates,
        matid_weak_seed_mantle
    )

    return matid
end

end # module LithosphericExtensionWeakSeed 