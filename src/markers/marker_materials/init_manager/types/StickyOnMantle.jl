module StickyOnMantle

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
    thick_air = model.geometry.parameters.sticky_air_geometry.thick_air.value

    matid_domains_dict = model.materials.dicts.matid_domains
    matid_asthenosphere = matid_domains_dict["Asthenosphere"]
    matid_sticky_water = matid_domains_dict["Ocean"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    for imarker in 1:marknum
        y_marker = marker_y[imarker]
        matid = define_material(
            y_marker, matid_sticky_water, matid_asthenosphere, thick_air
        )
        marker_matid[imarker] = matid
    end

    return nothing
end

function define_material(
    y_marker::Float64,
    matid_sticky_water::Int16,
    matid_asthenosphere::Int16,
    thick_air::Float64
)::Int16
    if y_marker <= thick_air
        matid = matid_sticky_water
    else
        matid = matid_asthenosphere
    end
    return matid
end

end # module StickyOnMantle 