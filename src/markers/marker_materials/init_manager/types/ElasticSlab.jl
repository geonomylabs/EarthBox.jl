module ElasticSlab

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

    domains = model.materials.dicts.matid_domains
    matid_fluid_around_elastic_slab_a = domains["FluidAroundElasticSlabA"]
    matid_fluid_around_elastic_slab_b = domains["FluidAroundElasticSlabB"]
    matid_elastic_slab_a = domains["ElasticSlabA"]
    matid_elastic_slab_b = domains["ElasticSlabB"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = define_material(
            x_marker, y_marker, xsize, ysize,
            matid_fluid_around_elastic_slab_a,
            matid_fluid_around_elastic_slab_b,
            matid_elastic_slab_a,
            matid_elastic_slab_b
        )
        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_material(
    x_marker::Float64,
    y_marker::Float64,
    xsize::Float64,
    ysize::Float64,
    matid_fluid_around_elastic_slab_a::Int16,
    matid_fluid_around_elastic_slab_b::Int16,
    matid_elastic_slab_a::Int16,
    matid_elastic_slab_b::Int16
)::Int16
    matid = matid_fluid_around_elastic_slab_a

    # Checkerboard structure
    m1_check = calculate_checkerboard_flag_m1(x_marker, xsize)
    m2_check = calculate_checkerboard_flag_m2(y_marker, ysize)

    if m2_check == m1_check
        matid = matid_fluid_around_elastic_slab_b
    end

    # Slab
    if in_elastic_slab(y_marker, x_marker, ysize, xsize)
        if matid == matid_fluid_around_elastic_slab_a
            matid = matid_elastic_slab_a
        else
            matid = matid_elastic_slab_b
        end
    end
    return matid
end

function calculate_checkerboard_flag_m1(x_marker::Float64, xsize::Float64)::Float64
    arg = x_marker/(xsize/20.0)
    m1_check = Float64(ceil(Int, arg)) + 1.0

    arg = m1_check/2
    m1_check = m1_check - Float64(ceil(Int, arg))*2.0
    return m1_check
end

function calculate_checkerboard_flag_m2(y_marker::Float64, ysize::Float64)::Float64
    arg = y_marker/(ysize/20.0)
    m2_check = Float64(ceil(Int, arg)) + 1.0

    arg = m2_check/2
    m2_check = m2_check - Float64(ceil(Int, arg))*2.0
    return m2_check
end

function in_elastic_slab(
    y_marker::Float64,
    x_marker::Float64,
    ysize::Float64,
    xsize::Float64
)::Bool
    ymin_nondim, ymax_nondim, xmax_nondim = get_elastic_slab_non_dimensional_limits()
    check = false
    if ymin_nondim < y_marker/ysize < ymax_nondim && x_marker/xsize < xmax_nondim
        check = true
    end
    return check
end

function get_elastic_slab_non_dimensional_limits()::Tuple{Float64, Float64, Float64}
    ymin_nondim = 0.2
    ymax_nondim = 0.8
    xmax_nondim = 0.8
    return (ymin_nondim, ymax_nondim, xmax_nondim)
end

end # module ElasticSlab 