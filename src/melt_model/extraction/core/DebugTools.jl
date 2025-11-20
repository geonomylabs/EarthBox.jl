module DebugTools

import EarthBox.ModelDataContainer: ModelData

function initialize_debug_variable_arrays(nmarkers_magma_mantle::Int)
    yshallow_partial_melt_min = 1e32
    yshallow_partial_melt_max = -1e32

    yshallow_mantle_min = 1e32
    yshallow_mantle_max = -1e32

    mantle_injection_search_xmin_min = 1e32
    mantle_injection_search_xmin_max = -1e32

    mantle_injection_search_xmax_min = 1e32
    mantle_injection_search_xmax_max = -1e32

    sub_domain_index_min = 9999
    sub_domain_index_max = -9999

    debug_variables_floats = [
        yshallow_partial_melt_min,
        yshallow_partial_melt_max,
        yshallow_mantle_min,
        yshallow_mantle_max,
        mantle_injection_search_xmin_min,
        mantle_injection_search_xmin_max,
        mantle_injection_search_xmax_min,
        mantle_injection_search_xmax_max
    ]

    sub_domain_index = -1
    icount_sub_domain_index = 0

    debug_variables_ints = [
        sub_domain_index,
        icount_sub_domain_index,
        sub_domain_index_min,
        sub_domain_index_max
    ]

    sub_domain_index_array = fill(-1, nmarkers_magma_mantle)

    return debug_variables_floats, debug_variables_ints, sub_domain_index_array
end

function update_subdomain_index_array(
    iuse_shallow_mantle_injection::Int,
    debug_variables_ints::Vector{Int},
    sub_domain_index::Int,
    sub_domain_index_array::Vector{Int}
)
    if iuse_shallow_mantle_injection == 1
        icount_sub_domain_index = debug_variables_ints[2]
        sub_domain_index_array[icount_sub_domain_index] = sub_domain_index
        icount_sub_domain_index += 1
        debug_variables_ints[2] = icount_sub_domain_index
    end
end

function update_debugging_arrays(
    imarker_mantle_shallow::Int,
    yshallow_mantle::Float64,
    mantle_injection_search_xmin::Float64,
    mantle_injection_search_xmax::Float64,
    sub_domain_index::Int,
    debug_variables_floats::Vector{Float64},
    debug_variables_ints::Vector{Int}
)
    if imarker_mantle_shallow != -999
        debug_variables_floats[3], debug_variables_floats[4] = 
            update_min_max_value(yshallow_mantle, debug_variables_floats[3], 
                               debug_variables_floats[4])

        debug_variables_floats[5], debug_variables_floats[6] = 
            update_min_max_value(mantle_injection_search_xmin,
                               debug_variables_floats[5],
                               debug_variables_floats[6])

        debug_variables_floats[7], debug_variables_floats[8] = 
            update_min_max_value(mantle_injection_search_xmax,
                               debug_variables_floats[7],
                               debug_variables_floats[8])

        debug_variables_ints[3], debug_variables_ints[4] = 
            update_min_max_value_int(sub_domain_index,
                                   debug_variables_ints[3],
                                   debug_variables_ints[4])
    end
end

function check_marker_y_coordinate(
    model::ModelData,
    marker_index::Int,
    yshallow_max::Float64;
    y_limit::Float64=25_000.0
)
    x = model.markers.arrays.location.marker_x.array[marker_index]
    y = model.markers.arrays.location.marker_y.array[marker_index]
    
    if y >= y_limit
        println("!!! Warning !!!: Marker y-coordinate ", y, " x-coordinate ", x,
                " with index ", marker_index,
                " exceeds limit y_limit ", y_limit)
    end
    
    if y > yshallow_max
        println("!!! Warning !!!: Marker y-coordinate ", y, " x-coordinate ", x,
                " with index ", marker_index,
                " exceeds maximum emplaced y-coordinate ", yshallow_max,
                " This does not make sense")
    end
end

function update_min_max_value(
    value::Float64,
    value_min::Float64,
    value_max::Float64
)
    value_min = min(value_min, value)
    value_max = max(value_max, value)
    return value_min, value_max
end

function update_min_max_value_int(
    value::Int,
    value_min::Int,
    value_max::Int
)
    value_min = min(value_min, value)
    value_max = max(value_max, value)
    return value_min, value_max
end

function print_debug_for_magma_body(
    magma_height_limit::Float64,
    iuse_shallow_mantle_injection::Int,
    sub_domain_index_array::Vector{Int},
    debug_variables_floats::Vector{Float64},
    debug_variables_ints::Vector{Int}
)
    print_yshallow_min_max(
        debug_variables_floats[1],
        debug_variables_floats[2],
        "partial melt replacement"
    )
    
    if iuse_shallow_mantle_injection == 1
        filtered_array = filter(x -> x != -1, sub_domain_index_array)
        if !isempty(filtered_array)
            sub_domain_index_mode = mode(filtered_array)
            sub_domain_index_mean = mean(filtered_array)
            sub_domain_index_std = std(filtered_array)
            print_sub_domain_stats(
                sub_domain_index_mode,
                sub_domain_index_mean,
                sub_domain_index_std
            )
        end
        
        print_mantle_injection_x_min_max_info(
            debug_variables_floats[5],
            debug_variables_floats[6],
            debug_variables_floats[7],
            debug_variables_floats[8],
            debug_variables_ints[3],
            debug_variables_ints[4]
        )
        
        print_yshallow_min_max(
            debug_variables_floats[3],
            debug_variables_floats[4],
            "mantle replacement"
        )
        
        print_emplacement_height(
            debug_variables_floats[3],
            debug_variables_floats[4],
            magma_height_limit
        )
        
        check_emplacement_height(
            debug_variables_floats[3],
            debug_variables_floats[4],
            magma_height_limit
        )
    end
end

function print_yshallow_min_max(
    yshallow_min::Float64,
    yshallow_max::Float64,
    description::String
)
    println(">> Minimum y-coordinate for ", description, " ", yshallow_min)
    println(">> Maximum y-coordinate for ", description, " ", yshallow_max)
end

function print_sub_domain_stats(
    sub_domain_index_mode::Float64,
    sub_domain_index_mean::Float64,
    sub_domain_index_std::Float64
)
    println(">> Subdomain index mode ", sub_domain_index_mode)
    println(">> Subdomain index mean ", sub_domain_index_mean)
    println(">> Subdomain index standard deviation ", sub_domain_index_std)
end

function print_mantle_injection_x_min_max_info(
    injection_domain_xmin_min::Float64,
    injection_domain_xmin_max::Float64,
    injection_domain_xmax_min::Float64,
    injection_domain_xmax_max::Float64,
    sub_domain_index_min::Int,
    sub_domain_index_max::Int
)
    println(">> Minimum left x-coordinate for injection domain ", 
            injection_domain_xmin_min)
    println(">> Maximum left x-coordinate for injection domain ", 
            injection_domain_xmin_max)
    println(">> Minimum right x-coordinate for injection domain ", 
            injection_domain_xmax_min)
    println(">> Maximum right x-coordinate for injection domain ", 
            injection_domain_xmax_max)
    println(">> Minimum subdomain index ", sub_domain_index_min)
    println(">> Maximum subdomain index ", sub_domain_index_max)
end

function print_emplacement_height(
    yshallow_mantle_min::Float64,
    yshallow_mantle_max::Float64,
    magma_height_limit::Float64
)
    emplacement_height = yshallow_mantle_max - yshallow_mantle_min
    println(">> Current emplacement height ", emplacement_height,
            " vs. magma height limit ", magma_height_limit)
    println("")
end

function check_emplacement_height(
    yshallow_mantle_min::Float64,
    yshallow_mantle_max::Float64,
    magma_height_limit::Float64
)
    emplacement_height = yshallow_mantle_max - yshallow_mantle_min
    if emplacement_height > 2.0 * magma_height_limit
        println("!!! Warning !!!: Emplacement height ", emplacement_height,
                " exceeds 2 x magma height limit of ", magma_height_limit)
        println("")
    end
end

function check_emplacement_height_in_loop(
    model::ModelData,
    imarker_shallow::Int,
    yshallow_mantle::Float64,
    yshallow_mantle_min::Float64,
    mantle_injection_search_xmin::Float64,
    mantle_injection_search_xmax::Float64,
    xshallow_partial_melt::Float64,
    injection_width::Float64;
    threshold::Float64=8000.0
)
    x = model.markers.arrays.location.marker_x.array[imarker_shallow]
    y = model.markers.arrays.location.marker_y.array[imarker_shallow]

    emplacement_height = yshallow_mantle - yshallow_mantle_min
    if emplacement_height > threshold
        println("!!! Warning !!!: Emplacement height is greater than threshold: ",
                emplacement_height, " exceeds threshold of ", threshold)
        println(">> mantle_injection_search_xmin ", mantle_injection_search_xmin,
                " mantle_injection_search_xmax ", mantle_injection_search_xmax)
        println(">> xshallow_partial_melt ", xshallow_partial_melt,
                " injection_width ", injection_width)
        println(">> imarker ", imarker_shallow, " x ", x, " y ", y,
                " yshallow_mantle ", yshallow_mantle)
        println("")
    end
end

end # module 