module SimpleSediment

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
    marker_y = model.markers.arrays.location.marker_y.array

    thick_air = model.geometry.parameters.sticky_air_geometry.thick_air.value
    thick_sediment_meters = 5_000.0

    domains = model.materials.dicts.matid_domains
    matid_sticky_air = domains["Atmosphere"]
    matid_sediment = domains["SedimentaryBasin"]
    matid_asthenosphere = domains["Asthenosphere"]

    marknum = model.markers.parameters.distribution.marknum.value
    for imarker in 1:marknum
        y_marker = marker_y[imarker]
        if y_marker <= thick_air
            matid = matid_sticky_air
        elseif thick_air < y_marker <= thick_air + thick_sediment_meters
            matid = matid_sediment
        else
            matid = matid_asthenosphere
        end
        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

end # module SimpleSediment 