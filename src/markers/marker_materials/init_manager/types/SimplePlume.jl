module SimplePlume

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import ..Plume: in_plume_region_check
import ..InitManager.InitStructs: PlumeGeometry

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
    thick_air = model.geometry.parameters.sticky_air_geometry.thick_air.value

    plume = model.geometry.parameters.plume
    plume_geometry = PlumeGeometry(
        plume_radius = plume.plume_radius.value,
        plume_center_x = plume.plume_center_x.value,
        plume_center_y = plume.plume_center_y.value,
        plume_head_thick = plume.plume_head_thick.value
    )

    matid_domains_dict = model.materials.dicts.matid_domains
    matid_sticky_air = matid_domains_dict["Atmosphere"]
    matid_asthenosphere = matid_domains_dict["Asthenosphere"]
    matid_plume = matid_domains_dict["MantlePlume"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = define_material(
            x_marker,
            y_marker,
            thick_air,
            matid_sticky_air,
            matid_asthenosphere,
            matid_plume,
            plume_geometry
        )
        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_material(
    x_marker::Float64,
    y_marker::Float64,
    thick_air::Float64,
    matid_sticky_air::Int16,
    matid_asthenosphere::Int16,
    matid_plume::Int16,
    plume_geometry::PlumeGeometry
)::Int16
    matid = matid_asthenosphere
    if y_marker < thick_air
        matid = matid_sticky_air
    end
    if in_plume_region_check(x_marker, y_marker, plume_geometry)
        matid = matid_plume
    end
    return matid
end

end # module SimplePlume 