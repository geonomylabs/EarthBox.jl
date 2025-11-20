module FindShallowest

import EarthBox.MathTools: generate_normal_random_number

""" Calculate index and depth of the shallowest partially molten marker.

# Returns
- imarker_shallow::Int
    - Index of shallowest partially molten marker.
- yshallow::Float64
    - Y-coordinate of shallowest partially molten marker.
"""
function find_shallowest_partially_molten_mantle_marker_opt(
    marker_matid::Vector{Int16},
    marker_y::Vector{Float64},
    mantle_melting_mat_ids::Vector{Int16},
    layer_counts::Vector{Int64},
    marker_indices::Vector{Vector{Int64}}
)::Tuple{Int, Float64}
    imarker_shallow = -999
    yshallow = 1e32
    found_marker = false
    nlayers = length(layer_counts)
    # Consider optimizing this for column-major 2D arrays
    for ilayer in 1:nlayers
        nmarkers_layer = layer_counts[ilayer]
        for j in 1:nmarkers_layer
            # get index of partially molten mantle marker
            imarker = marker_indices[ilayer][j]
            matid = marker_matid[imarker]
            # Only consider ID's associated with the mantle melting model
            if matid in mantle_melting_mat_ids
                found_marker = true
                # Get coordinates and material ID of marker
                y_marker = marker_y[imarker]
                if y_marker < yshallow
                    yshallow = y_marker
                    imarker_shallow = imarker
                end
            end
        end
        if found_marker
            break
        end
    end
    return imarker_shallow, yshallow
end

""" Calculate index and depth of the shallowest mantle marker in injection domain.

Find the shallowest marker in the injection domain with material id 
in mantle_mat_ids. The injection domain is a column of the model centered on
the shallowest partially molten mantle marker with a width equal to the
injection width. Only markers with indices in the marker mantle injection 
search domain array are searched for efficiency. The injection search domain
is has a width equal to the injection width times a user specified factor
(i.e. it is larger than the injection domain by for example a factor of 
three).

# Returns
- imarker_shallow::Int
    - Index of shallowest partially molten marker.
- yshallow::Float64
    - Y-coordinate of shallowest partially molten marker.
- xmin::Float64
    - Minimum x-coordinate of the subdomain.
- xmax::Float64
    - Maximum x-coordinate of the subdomain.
- sub_domain_index::Int
    - Index of the subdomain.
"""
function find_shallowest_mantle_marker_random_opt(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    marker_indices_mantle_injection_search_domain::Vector{Int64},
    mantle_mat_ids::Vector{Int16},
    x_shallow_partial_melt::Float64,
    injection_width::Float64;
    iuse_random_injection_subdomain::Int=0,
    iuse_normal_injection_subdomain::Int=0,
    number_of_injection_subdomains::Int=10
)::Tuple{Int, Float64, Float64, Float64, Int}
    xmin = x_shallow_partial_melt - injection_width/2.0
    xmax = x_shallow_partial_melt + injection_width/2.0
    sub_domain_index = -1

    if iuse_random_injection_subdomain == 1 || iuse_normal_injection_subdomain == 1
        xmin, xmax, sub_domain_index = calculate_random_subdomain_limits(
            injection_width, xmin,
            number_of_injection_subdomains,
            iuse_normal_injection_subdomain
        )
    end

    imarker_shallow, yshallow = find_shallowest_marker_in_mantle_injection_domain(
        marker_x,
        marker_y,
        marker_matid,
        marker_indices_mantle_injection_search_domain,
        mantle_mat_ids,
        xmin,
        xmax
    )

    return imarker_shallow, yshallow, xmin, xmax, sub_domain_index
end

""" Calculate random subdomain limits.

# Returns
- xmin_subdomain::Float64
    - Minimum x-coordinate of the subdomain.
- xmax_subdomain::Float64
    - Maximum x-coordinate of the subdomain.
- random_sub_domain_index::Int
    - Index of the random subdomain.
"""
function calculate_random_subdomain_limits(
    injection_width::Float64,
    xmin_main_domain::Float64,
    number_of_injection_subdomains::Int=10,
    iuse_normal_injection_subdomain::Int=0
)::Tuple{Float64, Float64, Int}
    dx = injection_width/number_of_injection_subdomains
    if iuse_normal_injection_subdomain == 1
        random_sub_domain_index = calculate_normal_injection_subdomain(
            0, number_of_injection_subdomains)
    else
        random_sub_domain_index = rand(0:number_of_injection_subdomains-1)
    end
    xmin_subdomain = xmin_main_domain + random_sub_domain_index * dx
    xmax_subdomain = xmin_subdomain + dx
    return xmin_subdomain, xmax_subdomain, random_sub_domain_index
end

""" Calculate the injection location using a normal distribution. """
function calculate_normal_injection_subdomain(
    minimum_domain_number::Int64,
    maximum_domain_number::Int64
)::Int64
    mean = minimum_domain_number + maximum_domain_number / 2.0
    std_dev = maximum_domain_number / 4.0
    subdomain_number = generate_normal_random_number(
        mean=mean, standard_deviation=std_dev)
    while subdomain_number < minimum_domain_number || 
          subdomain_number > minimum_domain_number + maximum_domain_number
        subdomain_number = generate_normal_random_number(
            mean=mean, standard_deviation=std_dev)
    end
    return floor(Int, subdomain_number)
end

""" Calculate index and depth of the shallowest marker.

# Returns
- imarker_shallow::Int
    - Index of shallowest marker.
- yshallow::Float64
    - Y-coordinate of shallowest marker.
"""
function find_shallowest_marker_in_mantle_injection_domain(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    marker_indices_mantle_injection_domain::Vector{Int64},
    mat_ids::Vector{Int16},
    xmin::Float64,
    xmax::Float64
)::Tuple{Int, Float64}
    nmarkers_injection_domain = length(marker_indices_mantle_injection_domain)
    imarker_shallow = -999
    yshallow = 1e32
    for i in 1:nmarkers_injection_domain
        imarker = marker_indices_mantle_injection_domain[i]
        x_marker = marker_x[imarker]
        if xmin <= x_marker <= xmax
            matid = marker_matid[imarker]
            if matid in mat_ids
                y_marker = marker_y[imarker]
                if y_marker < yshallow
                    yshallow = y_marker
                    imarker_shallow = imarker
                end
            end
        end
    end
    return imarker_shallow, yshallow
end

""" Calculate index and depth of the shallowest marker.

# Returns
- imarker_shallow::Int
    - Index of shallowest marker.
- yshallow::Float64
    - Y-coordinate of shallowest marker.
"""
function find_shallowest_marker(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    mat_ids::Vector{Int16},
    xmin::Float64,
    xmax::Float64,
    ymax::Float64
)::Tuple{Int, Float64}
    nmarkers = length(marker_y)
    imarker_shallow = -999
    yshallow = 1e32
    for imarker in 1:nmarkers
        x_marker = marker_x[imarker]
        if xmin <= x_marker <= xmax
            y_marker = marker_y[imarker]
            if y_marker < ymax
                matid = marker_matid[imarker]
                if matid in mat_ids
                    if y_marker < yshallow
                        yshallow = y_marker
                        imarker_shallow = imarker
                    end
                end
            end
        end
    end
    return imarker_shallow, yshallow
end

""" Calculate index and depth of the shallowest partially molten marker.

# Returns
- imarker_shallow::Int
    - Index of shallowest partially molten marker.
- yshallow::Float64
    - Y-coordinate of shallowest partially molten marker.
"""
function find_shallowest_partially_molten_mantle_marker(
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    nmarkers_partial_melt::Int,
    mantle_melting_mat_ids::Vector{Int16},
    partial_melt_flags::Vector{Int64}
)::Tuple{Int, Float64}
    imarker_shallow = -999
    yshallow = 1e32
    for j in 1:nmarkers_partial_melt
        # get index of partially molten mantle marker
        imarker = partial_melt_flags[j]
        # Get coordinates and material ID of marker
        y_marker = marker_y[imarker]
        matid = marker_matid[imarker]
        # Only consider ID's associated with the mantle melting model
        if matid in mantle_melting_mat_ids
            if y_marker < yshallow
                yshallow = y_marker
                imarker_shallow = imarker
            end
        end
    end
    return imarker_shallow, yshallow
end

""" Calculate index and coordinates shallowest magma marker.

Here magma refers to extracted gabbroic magma and partially molten gabbro.

# Returns
- imarker_shallow::Int
    - Index of shallowest partially molten marker.
- xshallow::Float64
    - X-coordinate of shallowest partially molten marker.
- yshallow::Float64
    - Y-coordinate of shallowest partially molten marker.
"""
function find_shallowest_gabbroic_partially_molten_or_magma_marker(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    molten_gabbro_ids::Vector{Int16};
    xstart::Float64=-1e39,
    xend::Float64=1e39
)::Tuple{Int, Float64, Float64}
    nmarkers = length(marker_x)
    imarker_shallow = -999
    yshallow = 1e32
    xshallow = 1e32
    for imarker in 1:nmarkers
        @inbounds x_marker = marker_x[imarker]
        if xstart <= x_marker <= xend
            @inbounds begin
                y_marker = marker_y[imarker]
                matid = marker_matid[imarker]
            end
            if matid in molten_gabbro_ids
                if y_marker < yshallow
                    xshallow = x_marker
                    yshallow = y_marker
                    imarker_shallow = imarker
                end
            end
        end
    end
    return imarker_shallow, xshallow, yshallow
end

end # module 