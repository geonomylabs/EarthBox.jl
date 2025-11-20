module SeafloorSpreading

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!

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

    thick_upper_lith = earth_layering.thick_upper_lith.value
    thick_middle_lith = earth_layering.thick_middle_lith.value
    thick_lower_lith = earth_layering.thick_lower_lith.value
    thick_lithosphere = thick_upper_lith + thick_middle_lith + thick_lower_lith

    thick_air = model.geometry.parameters.sticky_air_geometry.thick_air.value

    matid_domains_dict = model.materials.dicts.matid_domains

    matid_sticky_air = matid_domains_dict["Atmosphere"]
    matid_lithosphere = matid_domains_dict["MantleLithosphere"]
    matid_asthenosphere = matid_domains_dict["Asthenosphere"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        y_marker = marker_y[imarker]
        matid = define_material(
            y_marker,
            thick_air,
            thick_lithosphere,
            matid_sticky_air,
            matid_lithosphere,
            matid_asthenosphere
        )
        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_material(
    y_marker::Float64,
    thick_air::Float64,
    thick_lithosphere::Float64,
    matid_sticky_air::Int16,
    matid_lithosphere::Int16,
    matid_asthenosphere::Int16
)::Int16
    if y_marker < thick_air
        matid = matid_sticky_air
    elseif thick_air <= y_marker <= thick_air + thick_lithosphere
        matid = matid_lithosphere
    else
        matid = matid_asthenosphere
    end
    return matid
end

end # module SeafloorSpreading 