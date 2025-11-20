module ViscousBlock

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
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = define_material(x_marker, y_marker, xsize, ysize)
        marker_matid[imarker] = matid
    end

    return nothing
end

function define_material(
    x_marker::Float64,
    y_marker::Float64,
    xsize::Float64,
    ysize::Float64
)::Int
    matid_fluid_around_block = 2
    matid_viscous_block = 3
    matid = matid_fluid_around_block
    dx_frac = x_marker/xsize
    dy_frac = y_marker/ysize
    if 0.4 <= dx_frac <= 0.6 && 0.1 <= dy_frac <= 0.3
        matid = matid_viscous_block
    end
    return matid
end

end # module ViscousBlock 