module Uniform

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
    domains = model.materials.dicts.matid_domains
    matid = domains["GeneralDomain"]
    marknum = model.markers.parameters.distribution.marknum.value
    marker_matid = model.markers.arrays.material.marker_matid.array

    for imarker in 1:marknum
        marker_matid[imarker] = Int16(matid)
    end

    return nothing
end

end # module Uniform 