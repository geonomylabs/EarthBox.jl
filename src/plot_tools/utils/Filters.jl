module Filters

import Printf: @printf
import EarthBox.ModelDataContainer: ModelData

function apply_dimensional_filter(
    marker_x_km::Vector{Float64},
    marker_y_km::Vector{Float64},
    marker_scalar::Union{Vector{Float32}, Vector{Float64}, Vector{Int16}, Vector{Int32}, Vector{Int64}},
    plot_dimensions::Tuple{Float64, Float64, Float64, Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    marker_x_km, marker_y_km, marker_scalar = filter_markers_for_plot_dimensions(
        marker_x_km, marker_y_km, marker_scalar, plot_dimensions)
    return marker_x_km, marker_y_km, marker_scalar
end

function filter_markers_for_plot_dimensions(
    marker_x_km::Vector{Float64},
    marker_y_km::Vector{Float64},
    marker_scalar::Union{Vector{Float32}, Vector{Float64}, Vector{Int16}, Vector{Int32}, Vector{Int64}},
    plot_dimensions::Tuple{Float64, Float64, Float64, Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    @printf("   >> Applying dimensional filter for plot dimensions: %f, %f, %f, %f\n", plot_dimensions...)
    nmarker = length(marker_x_km)
    @printf("      >> Number of markers before dimensional filtering: %d\n", nmarker)
    (
        marker_x_km_filtered_tmp,
        marker_y_km_filtered_tmp,
        marker_scalar_filtered_tmp
    ) = create_temporary_marker_arrays_simple(nmarker)

    nfilter = 0
    for i in 1:nmarker
        x = marker_x_km[i]
        y = marker_y_km[i]
        if (
            plot_dimensions[1] < x < plot_dimensions[2]
            && plot_dimensions[3] < y < plot_dimensions[4]
        )
            marker_x_km_filtered_tmp[nfilter + 1] = marker_x_km[i]
            marker_y_km_filtered_tmp[nfilter + 1] = marker_y_km[i]
            marker_scalar_filtered_tmp[nfilter + 1] = marker_scalar[i]
            nfilter += 1
        end
    end

    (
        marker_x_km_filtered, marker_y_km_filtered,
        marker_scalar_filtered
    ) = clean_filtered_marker_arrays_simple(
        nfilter, marker_x_km_filtered_tmp,
        marker_y_km_filtered_tmp, marker_scalar_filtered_tmp
    )
    nmarker = length(marker_x_km_filtered)
    println("      >> Number of markers after dimensional filtering: $nmarker")

    return (
        marker_x_km_filtered, marker_y_km_filtered,
        marker_scalar_filtered
    )
end

function filter_markers_for_age_plot(
    marker_x_km::Vector{Float64},
    marker_y_km::Vector{Float64},
    marker_age::Vector{Float64},
    marker_matid::Union{Vector{Float64}, Vector{Int16}},
    matids_to_keep::Vector{Int16}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int16}}
    nmarker = length(marker_x_km)
    
    (
        marker_x_km_filtered_tmp,
        marker_y_km_filtered_tmp,
        marker_age_filtered_tmp,
        marker_matid_filtered_tmp
    ) = create_temporary_marker_arrays(nmarker)

    nfilter = 0
    for i in 1:nmarker
        mitype = marker_matid[i]
        if mitype in matids_to_keep
            marker_age_filtered_tmp[nfilter + 1] = marker_age[i]
            marker_x_km_filtered_tmp[nfilter + 1] = marker_x_km[i]
            marker_y_km_filtered_tmp[nfilter + 1] = marker_y_km[i]
            marker_matid_filtered_tmp[nfilter + 1] = mitype
            nfilter += 1
        end
    end

    (
        marker_x_km_filtered, marker_y_km_filtered,
        marker_age_filtered, marker_matid_filtered
    ) = clean_filtered_marker_arrays(
        nfilter, marker_age_filtered_tmp, marker_x_km_filtered_tmp,
        marker_y_km_filtered_tmp, marker_matid_filtered_tmp
    )

    return (
        marker_x_km_filtered, marker_y_km_filtered,
        marker_age_filtered, marker_matid_filtered
    )
end

function filter_markers_based_on_minimum_value(
    minimum_value::Float64,
    marker_x_km::Vector{Float64},
    marker_y_km::Vector{Float64},
    marker_scalar_array::Union{Vector{Float64}, Vector{Int16}},
    marker_matid::Union{Vector{Float64}, Vector{Int16}},
    matids_to_keep::Vector{Int16}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int16}}
    nmarker = length(marker_x_km)
    (
        marker_x_km_filtered_tmp,
        marker_y_km_filtered_tmp,
        marker_scalar_array_filtered_tmp,
        marker_matid_filtered_tmp
    ) = create_temporary_marker_arrays(nmarker)
    
    nfilter = 0
    for i in 1:nmarker
        matid = Int64(marker_matid[i])
        marker_scalar = marker_scalar_array[i]
        apply_filter = check_apply_minimum_filter(
            marker_scalar, minimum_value, matid, matids_to_keep)
        if apply_filter
            marker_scalar_array_filtered_tmp[nfilter + 1] = marker_scalar
            marker_x_km_filtered_tmp[nfilter + 1] = marker_x_km[i]
            marker_y_km_filtered_tmp[nfilter + 1] = marker_y_km[i]
            marker_matid_filtered_tmp[nfilter + 1] = matid
            nfilter += 1
        end
    end

    (
        marker_x_km_filtered,
        marker_y_km_filtered,
        marker_scalar_array_filtered,
        marker_matid_filtered
    ) = clean_filtered_marker_arrays(
        nfilter,
        marker_scalar_array_filtered_tmp,
        marker_x_km_filtered_tmp,
        marker_y_km_filtered_tmp,
        marker_matid_filtered_tmp
    )

    return (
        marker_x_km_filtered,
        marker_y_km_filtered,
        marker_scalar_array_filtered,
        marker_matid_filtered
    )
end

function check_apply_minimum_filter(
    marker_scalar::Float64,
    minimum_value::Float64,
    matid::Int64,
    matids_to_keep::Vector{Int16}
)::Bool
    check = false
    if marker_scalar > minimum_value
        check = true
    end
    
    max_matid = maximum(matids_to_keep)
    if max_matid > 0
        if !(matid in matids_to_keep)
            check = false
        end
    end

    return check
end

function strain_condition(
    marker_scalar::Float64,
    minimum_value::Float64,
    mitype::Int64,
    ids_to_exclude::Tuple{Int64, Int64},
    scalar_gt_zero::Float64
)::Bool
    check = false
    if (
        marker_scalar > minimum_value
        && !(mitype in ids_to_exclude)
        && scalar_gt_zero > 0
    )
        check = true
    end
    return check
end

function filter_markers_for_strain_plot(
    plastic_strain_minimum::Float64,
    marker_x_km::Vector{Float64},
    marker_y_km::Vector{Float64},
    marker_pfailure::Vector{Float64},
    marker_strain_plastic::Vector{Float64},
    marker_matid::Union{Vector{Float64}, Vector{Int16}},
    sticky_matids::Tuple{Int16, Int16},
    sed_matid::Int16,
    basalt_matid::Int16
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int16}}
    nmarker = length(marker_x_km)
    (
        marker_x_km_filtered_tmp,
        marker_y_km_filtered_tmp,
        marker_strain_plastic_filtered_tmp,
        marker_matid_filtered_tmp
    ) = create_temporary_marker_arrays(nmarker)

    nfilter = 0
    for i in 1:nmarker
        mitype = marker_matid[i]
        strain = marker_strain_plastic[i]
        plastic_failure_flag = marker_pfailure[i]
        if (
            strain > plastic_strain_minimum
            && !(mitype in sticky_matids)
            && mitype != sed_matid
            && mitype != basalt_matid
            && plastic_failure_flag > 0
        )
            marker_strain_plastic_filtered_tmp[nfilter + 1] = strain
            marker_x_km_filtered_tmp[nfilter + 1] = marker_x_km[i]
            marker_y_km_filtered_tmp[nfilter + 1] = marker_y_km[i]
            marker_matid_filtered_tmp[nfilter + 1] = mitype
            nfilter += 1
        end
    end

    (
        marker_x_km_filtered,
        marker_y_km_filtered,
        marker_strain_plastic_filtered,
        marker_matid_filtered
    ) = clean_filtered_marker_arrays(
        nfilter,
        marker_strain_plastic_filtered_tmp,
        marker_x_km_filtered_tmp,
        marker_y_km_filtered_tmp,
        marker_matid_filtered_tmp
    )

    return (
        marker_x_km_filtered,
        marker_y_km_filtered,
        marker_strain_plastic_filtered,
        marker_matid_filtered
    )
end

function create_temporary_marker_arrays(nmarker::Int64)
    marker_x_tmp = Vector{Float64}(undef, nmarker)
    # zeros(Float64, nmarker)
    marker_y_tmp = Vector{Float64}(undef, nmarker)
    # zeros(Float64, nmarker)
    marker_scalar_tmp = Vector{Float64}(undef, nmarker)
    # zeros(Float64, nmarker)
    marker_matid_tmp = Vector{Int16}(undef, nmarker)
    # zeros(Int64, nmarker)
    return marker_x_tmp, marker_y_tmp, marker_scalar_tmp, marker_matid_tmp
end

function clean_filtered_marker_arrays(
    nfilter::Int64,
    marker_scalar_tmp::Vector{Float64},
    marker_x_tmp::Vector{Float64},
    marker_y_tmp::Vector{Float64},
    marker_matid_tmp::Vector{Int16}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int16}}
    marker_scalar_filtered = Vector{Float64}(undef, nfilter)
    # zeros(Float64, nfilter)
    marker_x_filtered = Vector{Float64}(undef, nfilter)
    # zeros(Float64, nfilter)
    marker_y_filtered = Vector{Float64}(undef, nfilter)
    # zeros(Float64, nfilter)
    marker_matid_filtered = Vector{Int16}(undef, nfilter)
    # zeros(Int64, nfilter)
    
    for i in 1:nfilter
        marker_scalar_filtered[i] = marker_scalar_tmp[i]
        marker_x_filtered[i] = marker_x_tmp[i]
        marker_y_filtered[i] = marker_y_tmp[i]
        marker_matid_filtered[i] = marker_matid_tmp[i]
    end
    
    return (
        marker_x_filtered,
        marker_y_filtered,
        marker_scalar_filtered,
        marker_matid_filtered
    )
end

function create_temporary_marker_arrays_simple(
    nmarker::Int64
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    marker_x_tmp = Vector{Float64}(undef, nmarker)
    #zeros(Float64, nmarker)
    marker_y_tmp = Vector{Float64}(undef, nmarker)
    # zeros(Float64, nmarker)
    marker_scalar_tmp = Vector{Float64}(undef, nmarker)
    # zeros(Float64, nmarker)
    return marker_x_tmp, marker_y_tmp, marker_scalar_tmp
end

function clean_filtered_marker_arrays_simple(
    nfilter::Int64,
    marker_x_tmp::Vector{Float64},
    marker_y_tmp::Vector{Float64},
    marker_scalar_tmp::Union{Vector{Float64}, Vector{Int16}}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    marker_x_filtered = Vector{Float64}(undef, nfilter)
    # zeros(Float64, nfilter)
    marker_y_filtered = Vector{Float64}(undef, nfilter)
    # zeros(Float64, nfilter)
    marker_scalar_filtered = Vector{Float64}(undef, nfilter)
    # zeros(Float64, nfilter)
    
    for i in 1:nfilter
        marker_x_filtered[i] = marker_x_tmp[i]
        marker_y_filtered[i] = marker_y_tmp[i]
        marker_scalar_filtered[i] = marker_scalar_tmp[i]
    end
    
    return (
        marker_x_filtered,
        marker_y_filtered,
        marker_scalar_filtered
    )
end

end # module
