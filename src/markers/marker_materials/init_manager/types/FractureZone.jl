module FractureZone

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!

struct FractureZoneMaterialIDs
    matid_sticky_air::Int16
    matid_sediment::Int16
    matid_basaltic_oceanic_crust::Int16
    matid_gabbroic_oceanic_crust::Int16
    matid_mantle_lithosphere::Int16
    matid_weak_mantle_lithosphere::Int16
    matid_asthenosphere::Int
end

struct FractureZoneGeometry
    sticky_air_thickness::Float64
    x_fracture_zone_start::Float64
    x_fracture_zone_end::Float64
    base_of_sediments::Float64
    base_of_basaltic_crust::Float64
    base_of_gabbroic_crust::Float64
    base_of_younger_lithosphere::Float64
    base_of_older_lithosphere::Float64
    base_of_weak_lithosphere::Float64
end

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
    material_ids = define_material_ids(model)
    geometry = define_geometry(model)

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = define_material(x_marker, y_marker, material_ids, geometry)
        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_material_ids(model::ModelData)::FractureZoneMaterialIDs
    domains = model.materials.dicts.matid_domains
    return FractureZoneMaterialIDs(
        domains["Atmosphere"],
        domains["SedimentaryBasin"],
        domains["BasalticOceanicCrust"],
        domains["GabbroicOceanicCrust"],
        domains["MantleLithosphere"],
        domains["WeakMantleLithosphere"],
        domains["Asthenosphere"]
    )
end

function define_geometry(model::ModelData)::FractureZoneGeometry
    sticky_air_geometry = model.geometry.parameters.sticky_air_geometry
    fracture_zone = model.geometry.parameters.fracture_zone

    sticky_air_thickness = sticky_air_geometry.thick_air.value
    sediment_thickness = fracture_zone.sediment_thickness.value
    basaltic_oceanic_crust_thickness = fracture_zone.basaltic_oceanic_crust_thickness.value
    gabbroic_oceanic_crust_thickness = fracture_zone.gabbroic_oceanic_crust_thickness.value
    thickness_of_younger_lithosphere = fracture_zone.thickness_of_younger_lithosphere.value
    thickness_of_older_lithosphere = fracture_zone.thickness_of_older_lithosphere.value
    thickness_of_weak_lithosphere = fracture_zone.thickness_of_weak_lithosphere.value
    x_fracture_zone_start = fracture_zone.x_fracture_zone_start.value
    x_fracture_zone_end = fracture_zone.x_fracture_zone_end.value

    base_of_sediments = sticky_air_thickness + sediment_thickness
    base_of_basaltic_crust = base_of_sediments + basaltic_oceanic_crust_thickness
    base_of_gabbroic_crust = base_of_basaltic_crust + gabbroic_oceanic_crust_thickness
    base_of_younger_lithosphere = base_of_gabbroic_crust + thickness_of_younger_lithosphere
    base_of_older_lithosphere = base_of_gabbroic_crust + thickness_of_older_lithosphere
    base_of_weak_lithosphere = base_of_gabbroic_crust + thickness_of_weak_lithosphere

    return FractureZoneGeometry(
        sticky_air_thickness,
        x_fracture_zone_start,
        x_fracture_zone_end,
        base_of_sediments,
        base_of_basaltic_crust,
        base_of_gabbroic_crust,
        base_of_younger_lithosphere,
        base_of_older_lithosphere,
        base_of_weak_lithosphere
    )
end

function define_material(
    x_marker::Float64,
    y_marker::Float64,
    material_ids::FractureZoneMaterialIDs,
    geometry::FractureZoneGeometry
)::Int
    matid = set_asthenosphere(material_ids)
    matid = set_sticky_air(matid, y_marker, geometry, material_ids)
    matid = set_sediments(matid, y_marker, geometry, material_ids)
    matid = set_basaltic_crust(matid, y_marker, geometry, material_ids)
    matid = set_gabbroic_crust(matid, y_marker, geometry, material_ids)
    matid = set_lithosphere(matid, x_marker, y_marker, geometry, material_ids)
    matid = set_weak_gabbroic_oceanic_crust_of_fracture_zone(
        matid, x_marker, y_marker, geometry, material_ids)
    matid = set_weak_mantle_of_fracture_zone(
        matid, x_marker, y_marker, geometry, material_ids)
    return matid
end

function set_asthenosphere(material_ids::FractureZoneMaterialIDs)::Int16
    return material_ids.matid_asthenosphere
end

function set_sticky_air(
    matid::Int16,
    y_marker::Float64,
    geometry::FractureZoneGeometry,
    material_ids::FractureZoneMaterialIDs
)::Int16
    if y_marker <= geometry.sticky_air_thickness
        matid = material_ids.matid_sticky_air
    end
    return matid
end

function set_sediments(
    matid::Int16,
    y_marker::Float64,
    geometry::FractureZoneGeometry,
    material_ids::FractureZoneMaterialIDs
)::Int16
    if geometry.sticky_air_thickness < y_marker <= geometry.base_of_sediments
        matid = material_ids.matid_sediment
    end
    return matid
end

function set_basaltic_crust(
    matid::Int16,
    y_marker::Float64,
    geometry::FractureZoneGeometry,
    material_ids::FractureZoneMaterialIDs
)::Int16
    if geometry.base_of_sediments < y_marker <= geometry.base_of_basaltic_crust
        matid = material_ids.matid_basaltic_oceanic_crust
    end
    return matid
end

function set_gabbroic_crust(
    matid::Int16,
    y_marker::Float64,
    geometry::FractureZoneGeometry,
    material_ids::FractureZoneMaterialIDs
)::Int16
    if geometry.base_of_basaltic_crust < y_marker <= geometry.base_of_gabbroic_crust
        matid = material_ids.matid_gabbroic_oceanic_crust
    end
    return matid
end

function set_lithosphere(
    matid::Int16,
    x_marker::Float64,
    y_marker::Float64,
    geometry::FractureZoneGeometry,
    material_ids::FractureZoneMaterialIDs
)::Int16
    check_younger_lithosphere = in_younger_lithosphere(
        x_marker, y_marker, geometry.x_fracture_zone_end,
        geometry.base_of_gabbroic_crust, geometry.base_of_younger_lithosphere
    )
    check_older_lithosphere = in_older_lithosphere(
        x_marker, y_marker, geometry.x_fracture_zone_end,
        geometry.base_of_gabbroic_crust, geometry.base_of_older_lithosphere
    )
    if check_younger_lithosphere || check_older_lithosphere
        matid = material_ids.matid_mantle_lithosphere
    end
    return matid
end

function set_weak_gabbroic_oceanic_crust_of_fracture_zone(
    matid::Int16,
    x_marker::Float64,
    y_marker::Float64,
    geometry::FractureZoneGeometry,
    material_ids::FractureZoneMaterialIDs
)::Int16
    check_weak_gabbroic_oceanic_crust = in_weak_gabbroic_oceanic_crust_of_fracture_zone(
        x_marker, y_marker,
        geometry.x_fracture_zone_start, geometry.x_fracture_zone_end,
        geometry.base_of_basaltic_crust, geometry.base_of_gabbroic_crust
    )
    if check_weak_gabbroic_oceanic_crust
        matid = material_ids.matid_basaltic_oceanic_crust
    end
    return matid
end

function set_weak_mantle_of_fracture_zone(
    matid::Int16,
    x_marker::Float64,
    y_marker::Float64,
    geometry::FractureZoneGeometry,
    material_ids::FractureZoneMaterialIDs
)::Int16
    check_weak_mantle = in_weak_mantle_of_fracture_zone(
        x_marker, y_marker,
        geometry.x_fracture_zone_start, geometry.x_fracture_zone_end,
        geometry.base_of_gabbroic_crust, geometry.base_of_weak_lithosphere
    )
    if check_weak_mantle
        matid = material_ids.matid_weak_mantle_lithosphere
    end
    return matid
end

function in_weak_gabbroic_oceanic_crust_of_fracture_zone(
    x_marker::Float64,
    y_marker::Float64,
    x_fracture_zone_start::Float64,
    x_fracture_zone_end::Float64,
    base_of_basaltic_crust::Float64,
    base_of_gabbroic_crust::Float64
)::Bool
    x_weak_oceanic_crust_end = get_x_end_of_weak_oceanic_crust(
        y_marker, x_fracture_zone_end, base_of_gabbroic_crust)

    check_gabbroic_oceanic_crust = in_gabbroic_oceanic_crust(
        y_marker, base_of_basaltic_crust, base_of_gabbroic_crust
    )

    return x_fracture_zone_start < x_marker < x_weak_oceanic_crust_end && check_gabbroic_oceanic_crust
end

function get_x_end_of_weak_oceanic_crust(
    y_marker::Float64,
    x_fracture_zone_end::Float64,
    base_of_gabbroic_crust::Float64
)::Float64
    x_extension_of_weak_crust = (base_of_gabbroic_crust - y_marker)*25.0
    return x_fracture_zone_end + x_extension_of_weak_crust
end

function in_gabbroic_oceanic_crust(
    y_marker::Float64,
    base_of_basaltic_crust::Float64,
    base_of_gabbroic_crust::Float64
)::Bool
    return base_of_basaltic_crust < y_marker <= base_of_gabbroic_crust
end

function in_younger_lithosphere(
    x_marker::Float64,
    y_marker::Float64,
    x_fracture_zone_end::Float64,
    base_of_gabbroic_crust::Float64,
    base_of_younger_lithosphere::Float64
)::Bool
    return x_marker < x_fracture_zone_end && 
           base_of_gabbroic_crust < y_marker <= base_of_younger_lithosphere
end

function in_older_lithosphere(
    x_marker::Float64,
    y_marker::Float64,
    x_fracture_zone_end::Float64,
    base_of_gabbroic_crust::Float64,
    base_of_older_lithosphere::Float64
)::Bool
    return x_marker >= x_fracture_zone_end && 
           base_of_gabbroic_crust < y_marker <= base_of_older_lithosphere
end

function in_weak_mantle_of_fracture_zone(
    x_marker::Float64,
    y_marker::Float64,
    x_fracture_zone_start::Float64,
    x_fracture_zone_end::Float64,
    base_of_gabbroic_crust::Float64,
    base_of_weak_lithosphere::Float64
)::Bool
    return x_fracture_zone_start < x_marker < x_fracture_zone_end && 
           base_of_gabbroic_crust < y_marker <= base_of_weak_lithosphere
end

end # module FractureZone 