module MagmaBody

import EarthBox.PrintFuncs: print_info
import EarthBox.PrintFuncs: print_warning
import EarthBox.ModelDataContainer: ModelData
import EarthBox.SurfaceProcesses: calculate_age_ma
import EarthBox.SedimentWaterInterface: get_depth
import EarthBox: FindShallowest
import ..PartiallyMoltenZone

""" Make magma body at the top of mantle partially molten zone or mantle domain.

                        mantle search domain  
                     <----------------------->
            __________________________________________
            crust
            __________________________________________
            mantle              mmm
                            <--------->
                             injection
                               domain 
                        
                               __o__
                             // xxx \\            
                            //       \\
                           //         \\ 
                          // partially \\
                       __//   molten    \\___
                        ^                 ^
                    drainage           drainage  
                     divide             divide

Figure 1. Schematic of the extraction and emplacement scheme used in this
function. Depending on iuse_use_shallow_mantle_injection, extracted melt is 
either emplaced at the shallowest point (denoted by o) in the partially 
molten mantle zone in the associated drainage basin or the shallowest mantle
marker in the injection domain. The injection domain is a column of the 
model centered on the shallowest partially molten mantle marker with a 
width equal to the injection width. The symbols 'x' denote magma emplaced 
at the shallowest partially molten mantle marker in the drainage basin. The
symbols 'm' denote magma emplaced at the shallowest mantle marker in the
injection domain. The mantle search domain defines the subset of the mantle
domain that is searched for the shallowest mantle marker in the injection
domain for improved computational efficiency. The mantle search width
should be larger than the injection width otherwise erroneous injection
may occur. A warning is printed if the injection domain is outside of the 
mantle search domain.

# Returns
- xshallow_partial_melt_avg::Float64
    - Average x-coordinate of the shallowest partially molten mantle marker in
    the current drainage basin. An average is calculated because this 
    x-coordinate changes as markers are converted to magma.
- yshallow_partial_melt_avg::Float64
    - Average y-coordinate of the shallowest partially molten mantle marker.
"""
function extract_partial_melt_and_make_magma_body(
    model::ModelData,
    mantle_emplacement_mat_ids::Vector{Int16},
    nmarkers_magma_mantle::Int,
    nmarkers_partial_melt::Int,
    partial_melt_marker_indices::Vector{Int64},
    injection_width::Float64
)::Tuple{Float64, Float64}
    mantle_search_width = model.melting.parameters.extraction.mantle_search_width.value
    number_of_injection_subdomains = model.melting.parameters.extraction.number_of_injection_subdomains.value
    iuse_random_injection_subdomain = model.melting.parameters.options.iuse_random_injection_subdomain.value
    iuse_normal_injection_subdomain = model.melting.parameters.options.iuse_normal_injection_subdomain.value
    iuse_shallow_mantle_injection = model.melting.parameters.options.iuse_shallow_mantle_injection.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    matid_types = model.materials.dicts.matid_types
    matid_magma = matid_types["ExtractedGabbroicMagma"][1]

    age_ma = calculate_age_ma(model)

    # Divide partially molten markers into layers for more efficient search
    layer_counts, layered_partially_molten_marker_indices = 
        PartiallyMoltenZone.construct_layered_partially_molten_arrays(
            marker_x, marker_y, marker_matid,
            nmarkers_partial_melt,
            mantle_emplacement_mat_ids,
            partial_melt_marker_indices
        )

    # Calculate a subset of the total marker domain that will be searched for
    # the shallowest emplacement marker for extracted magma.
    marker_indices_mantle_injection_search_domain, xshallow_partial_melt_domain_initial,
    yshallow_partial_melt_domain_initial, mantle_search_xmin, mantle_search_xmax = 
        calculate_marker_indices_mantle_search_domain(
            model,
            mantle_emplacement_mat_ids,
            layer_counts,
            layered_partially_molten_marker_indices,
            search_width=mantle_search_width
        )

    use_threshold_depth = false
    if use_threshold_depth
        # Increase injection width if top of partially molten zone is beyond a
        # threshold depth. This may no longer be needed now that a bug has been
        # removed whereby melt from drainage basins was being applied to each
        # drainage basin.
        injection_width = update_shallow_injection_width_using_submud_depth(
            model, xshallow_partial_melt_domain_initial,
            yshallow_partial_melt_domain_initial, injection_width
        )
    end

    t1 = time()
    xshallow_partial_melt_avg = 0.0
    yshallow_partial_melt_avg = 0.0
    for _i in 1:nmarkers_magma_mantle
        # Get index and depth of the shallowest partially molten mantle marker
        # for the current drainage basin. Since the shallowest partially molten
        # markers are converted to magma they will be skipped in the next
        # iteration (if iuse_shallow_mantle_injection = 0).
        imarker_shallow, yshallow = 
            FindShallowest.find_shallowest_partially_molten_mantle_marker_opt(
                marker_matid,
                marker_y,
                mantle_emplacement_mat_ids,
                layer_counts,
                layered_partially_molten_marker_indices
            )

        xshallow_partial_melt = marker_x[imarker_shallow]
        xshallow_partial_melt_avg += xshallow_partial_melt
        yshallow_partial_melt_avg += yshallow

        if iuse_shallow_mantle_injection == 1 && imarker_shallow != -999
            icount_iter = 0
            imarker_mantle_shallow = -999
            while imarker_mantle_shallow == -999 && icount_iter < number_of_injection_subdomains
                # Find the shallowest mantle marker in the injection domain. The
                # injection domain is a column of the model centered on the
                # shallowest partially molten mantle marker in the current
                # drainage basin with a width equal to the injection width. Only
                # markers with indices in the marker_mantle injection search
                # domain array are searched for efficiency. The search domain
                # has a width equal to the injection width times a user
                # specified factor (e.g. factor = 3). This eliminates the need
                # to search the entire mantle domain for the shallowest mantle
                # marker.
                (
                    imarker_mantle_shallow, 
                    _yshallow_mantle, 
                    mantle_injection_search_xmin,
                    mantle_injection_search_xmax, 
                    _sub_domain_index
                ) = FindShallowest.find_shallowest_mantle_marker_random_opt(
                        marker_x,
                        marker_y,
                        marker_matid,
                        marker_indices_mantle_injection_search_domain,
                        mantle_emplacement_mat_ids,
                        xshallow_partial_melt,
                        injection_width,
                        iuse_random_injection_subdomain=iuse_random_injection_subdomain,
                        iuse_normal_injection_subdomain=iuse_normal_injection_subdomain,
                        number_of_injection_subdomains=number_of_injection_subdomains
                    )

                check_search_domains(
                    mantle_search_xmin, mantle_search_xmax,
                    mantle_injection_search_xmin, mantle_injection_search_xmax
                )

                icount_iter += 1
            end
            imarker_shallow = imarker_mantle_shallow
        end

        # Transform the shallowest marker to magma (either partially molten or
        # mantle depending on iuse_shallow_mantle_injection).
        if imarker_shallow != -999
            transform_marker_to_magma(model, imarker_shallow, age_ma, matid_magma)
        end
    end

    t2 = time()
    print_info("Time taken to extract partial melt and make magma body: $(t2 - t1) seconds", level=2)

    if nmarkers_magma_mantle > 0
        xshallow_partial_melt_avg = xshallow_partial_melt_avg/nmarkers_magma_mantle
        yshallow_partial_melt_avg = yshallow_partial_melt_avg/nmarkers_magma_mantle
    else
        xshallow_partial_melt_avg = -1e39
        yshallow_partial_melt_avg = -1e39
    end

    return xshallow_partial_melt_avg, yshallow_partial_melt_avg
end

""" Calculate indices of markers in the mantle search domain.

The mantle search domain is defined as the region where extracted magma
could be emplaced. The domain is defined relative to the shallowest 
partially molten mantle marker. The domain is defined in the x-direction
and extends to the left and right of the shallowest partially molten mantle
but does not extend below a certain depth beneath the top of the partially
molten mantle zone. The domain includes mantle markers with material IDs
that are in mantle_emplacement_mat_ids that are within the lateral search
width relative to the shallowest partially molten mantle marker.

This domain is used to limit the search domain for magma emplacement
and improve computational efficiency.

# Returns
- marker_indices_mantle_injection_search_domain::Vector{Int64}
    - Indices of markers in the mantle injection domain.
- xshallow_partial_melt_domain::Float64
    - X-coordinate of the shallowest partially molten mantle marker.
- yshallow_partial_melt_domain::Float64
    - Y-coordinate of the shallowest partially molten mantle marker.
- xstart::Float64
    - Minimum x-coordinate of the search domain.
- xend::Float64
    - Maximum x-coordinate of the search domain.
"""
function calculate_marker_indices_mantle_search_domain(
    model::ModelData,
    mantle_emplacement_mat_ids::Vector{Int16},
    layer_counts::Vector{Int64},
    layered_partially_molten_marker_indices::Vector{Vector{Int64}};
    search_width::Float64=25_000.0
)::Tuple{Vector{Int64}, Float64, Float64, Float64, Float64}
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marknum = model.markers.parameters.distribution.marknum.value

    (
        imarker_shallow_partial_melt_domain, yshallow_partial_melt_domain
    ) = FindShallowest.find_shallowest_partially_molten_mantle_marker_opt(
            marker_matid,
            marker_y,
            mantle_emplacement_mat_ids,
            layer_counts,
            layered_partially_molten_marker_indices
        )
    if imarker_shallow_partial_melt_domain < 0
        print_warning(
            "imarker_shallow_partial_melt_domain < 0: $(imarker_shallow_partial_melt_domain)", 
            level=2
            )
    end
    # Start here: index of shallowest partially molten marker is negative
    # in this code and the python version. The python/Numba version doesn't
    # crash when this is inserted into the array. The Julia version does
    # for good reason. This need to be managed more effectively.

    if imarker_shallow_partial_melt_domain >= 1
        xshallow_partial_melt_domain = marker_x[imarker_shallow_partial_melt_domain]
    else
        xshallow_partial_melt_domain = -1e39
    end

    xstart = xshallow_partial_melt_domain - search_width/2.0
    xend = xshallow_partial_melt_domain + search_width/2.0

    # depth below top of partial melt domain to search for mantle markers
    y_buffer = 20_000.0
    ymax = yshallow_partial_melt_domain + y_buffer

    marker_indices_tmp = Vector{Int64}(undef, marknum) #fill(-1, marknum)
    Threads.@threads for imarker in 1:marknum
        marker_indices_tmp[imarker] = -1
    end

    icount = 0
    for imarker in 1:marknum
        matid = marker_matid[imarker]
        if matid in mantle_emplacement_mat_ids
            x_marker = marker_x[imarker]
            if xstart <= x_marker <= xend
                y_marker = marker_y[imarker]
                if y_marker < ymax
                    marker_indices_tmp[icount + 1] = imarker
                    icount += 1
                end
            end
        end
    end
    nmarkers_injection_domain = icount

    marker_indices_mantle_injection_search_domain = Vector{Int64}(undef, nmarkers_injection_domain) 
    Threads.@threads for i in 1:nmarkers_injection_domain
        marker_indices_mantle_injection_search_domain[i] = marker_indices_tmp[i]
    end

    return (
        marker_indices_mantle_injection_search_domain,
        xshallow_partial_melt_domain,
        yshallow_partial_melt_domain,
        xstart,
        xend
    )
end

""" Update shallow mantle injection flag using submud depth.

If the submud depth of the initial peak of the partially molten mantle zone
is greater than the maximum shallow injection depth, the shallow mantle
injection flag is set to 0.

# Returns
- iuse_shallow_mantle_injection::Int
    - Updated shallow mantle injection flag.
"""
function update_shallow_injection_using_submud_depth(
    model::ModelData,
    xshallow_partial_melt_domain_initial::Float64,
    yshallow_partial_melt_domain_initial::Float64,
    iuse_shallow_mantle_injection::Int
)::Int
    gridt = model.topography.arrays.gridt.array

    y_mudline = get_depth(xshallow_partial_melt_domain_initial, gridt)

    submud_depth_partial_melt_peak = (
        yshallow_partial_melt_domain_initial - y_mudline)

    maximum_shallow_injection_depth = 
        model.melting.parameters.extraction.maximum_shallow_injection_depth.value

    if submud_depth_partial_melt_peak > maximum_shallow_injection_depth
        iuse_shallow_mantle_injection = 0
    end
    return iuse_shallow_mantle_injection
end

""" Update shallow mantle injection width using submud depth.

If the submud depth of the local maximum of the partially molten mantle zone
is greater than the maximum shallow injection depth, the shallow mantle
injection flag is set to 0 and the injection width is increased by a factor
of 5.

# Returns
- injection_width::Float64
    - Updated injection width.
"""
function update_shallow_injection_width_using_submud_depth(
    model::ModelData,
    xshallow_partial_melt_domain_initial::Float64,
    yshallow_partial_melt_domain_initial::Float64,
    injection_width::Float64
)::Float64
    gridt = model.topography.arrays.gridt.array

    y_mudline = get_depth(xshallow_partial_melt_domain_initial, gridt)

    submud_depth_partial_melt_peak = (
        yshallow_partial_melt_domain_initial - y_mudline)

    maximum_shallow_injection_depth = 
        model.melting.parameters.extraction.maximum_shallow_injection_depth.value

    if submud_depth_partial_melt_peak > maximum_shallow_injection_depth
        injection_width = injection_width * 5.0
    end
    return injection_width
end

""" Check search domains for problems. """
function check_search_domains(
    mantle_search_xmin::Float64,
    mantle_search_xmax::Float64,
    mantle_injection_search_xmin::Float64,
    mantle_injection_search_xmax::Float64
)::Nothing
    # Check if injection domain is within the mantle search domain. If not,
    # print a warning.
    if mantle_injection_search_xmin < mantle_search_xmin || 
       mantle_injection_search_xmax > mantle_search_xmax
        println(
            ">> Warning: Injection domain is outside of the mantle search domain. ",
            "Increase mantle search domain width."
        )
        println(">> mantle_search_xmin ", mantle_search_xmin)
        println(">> mantle_search_xmax ", mantle_search_xmax)
        println(">> mantle_injection_search_xmin ", mantle_injection_search_xmin)
        println(">> mantle_injection_search_xmax ", mantle_injection_search_xmax)
    end
    return nothing
end

""" Transform marker to pure molten material and define event age.

This function transforms marker material id stored in array marker_matid
to the material id of molten mantle material (matid_magma). This function
also resets a variety of marker arrays due to the change in composition.

# Updated Arrays
- model.markers.arrays.stress
    - marker_sxx.array: Normal stress of markers in Pascals.
    - marker_sxy.array: Shear stress of markers in Pascals.
- model.markers.arrays.rheology
    - marker_eta.array: Effective viscosity of markers in Pa.s.
    - marker_fric.array: Friction coefficient (sine of friction angle) of markers.
    - marker_pfailure.array: Plastic failure flag of markers.
- model.markers.arrays.strain
    - marker_exx.array: Normal strain rate of markers in 1/s.
    - marker_exy.array: Shear strain rate of markers in 1/s.
    - marker_GII.array: Strain of markers.
    - marker_strain_plastic.array: Plastic strain of markers.
    - marker_sr_ratio.array: Ratio of marker strain rate calculated using stress 
      change and Maxwell model over strain rate interpolated form grid.
- model.materials.arrays
    - marker_matid.array: Material ID of markers.
    - marker_serpentinization.array: Serpentinization ratio of markers.
- model.markers.arrays.melt
    - marker_meltfrac.array: Melt fraction of markers.
    - marker_extracted_meltfrac.array: Extracted melt fraction of marker.
    - marker_extractable_meltfrac.array: Extractable melt fraction of marker.
- model.markers.arrays.strat
    - marker_age.array: Age of formation of markers in Ma.
- model.markers.arrays.thermal
    - marker_TK.array: Temperature of markers in Kelvin.
"""
function transform_marker_to_magma(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    matid_magma::Int16;
    emplacement_temperature_kelvin::Float64=1473.0
)::Nothing
    @inbounds begin
        # set to the ID of molten mantle material
        model.markers.arrays.material.marker_matid.array[imarker] = matid_magma
        # Reset strain rate ratio
        model.markers.arrays.strain.marker_sr_ratio.array[imarker] = 1
        # Reset strain
        model.markers.arrays.strain.marker_GII.array[imarker] = 0
        model.markers.arrays.strain.marker_strain_plastic.array[imarker] = 0
        # Reset stress
        model.markers.arrays.stress.marker_sxx.array[imarker] = 0.0
        model.markers.arrays.stress.marker_sxy.array[imarker] = 0.0
        # Reset viscosity
        model.markers.arrays.rheology.marker_eta.array[imarker] = 0.0
        # Reset strain rate
        model.markers.arrays.strain.marker_exx.array[imarker] = 0.0
        model.markers.arrays.strain.marker_exy.array[imarker] = 0.0
        # Melt fraction; start at 1 since this is magma
        model.markers.arrays.melt.marker_meltfrac.array[imarker] = 1.0
        # Set total extracted melt to zero
        model.markers.arrays.melt.marker_extracted_meltfrac.array[imarker] = 0.0
        # Reset incremental melt extraction
        model.markers.arrays.melt.marker_extractable_meltfrac.array[imarker] = 0.0
        # Reset friction coefficient
        model.markers.arrays.rheology.marker_fric.array[imarker] = 
            model.materials.arrays.mat_plastic.array[matid_magma, 3]
        # Define age of melting
        model.markers.arrays.strat.marker_age.array[imarker] = age_ma
        # Set serpentinization to zero
        model.markers.arrays.material.marker_serpentinization.array[imarker] = 0.0
        # Set temperature to emplacement temperature
        model.markers.arrays.thermal.marker_TK.array[imarker] = emplacement_temperature_kelvin
    end
    return nothing
end

end # module 