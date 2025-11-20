module RayleighTaylorInstability

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import ....RayleighTaylor

function initialize!(
    model::ModelData;
    parameters::Union{Dict{String, Any}, Nothing}=nothing,
    material_domain_ids::Union{Dict{String, Int}, Nothing}=nothing
)::Nothing
    set_model_data!(model, parameters, material_domain_ids)
    initialize_material_ids!(model)
    RayleighTaylor.perturb_marker_y!(model)
    return nothing
end

function initialize_material_ids!(model::ModelData)::Nothing
    depth_interface_h1 = model.geometry.parameters.rayleigh_taylor.depth_interface_h1.value

    domains = model.materials.dicts.matid_domains
    matid_layer1_rayleigh_taylor = domains["Layer1RayleighTaylor"]
    matid_layer2_rayleigh_taylor = domains["Layer2RayleighTaylor"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        y_marker = marker_y[imarker]
        matid = define_material(
            y_marker,
            depth_interface_h1,
            matid_layer1_rayleigh_taylor,
            matid_layer2_rayleigh_taylor
        )
        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_material(
    y_marker::Float64,
    depth_interface::Float64,
    matid_layer1_rayleigh_taylor::Int16,
    matid_layer2_rayleigh_taylor::Int16
)::Int16
    matid = matid_layer1_rayleigh_taylor
    if y_marker >= depth_interface
        matid = matid_layer2_rayleigh_taylor
    end
    return matid
end

end # module RayleighTaylorInstability 