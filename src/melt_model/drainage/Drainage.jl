module Drainage

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_layer_opt
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_layer_opt_sorted
import EarthBox.ModelStructureManager.TopAndBottom: calculate_search_radius
import EarthBox.ModelStructureManager.SmoothSurface: smooth_surface

function calculate_melt_drainage_divides!(
    model::ModelData
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    topo_gridx = model.topography.arrays.gridt.array[1, :]
    partial_melt_gridy = calculate_top_of_mantle_partial_melt_domain(model, topo_gridx)
    divides_x = calculate_drainage_divides(topo_gridx, partial_melt_gridy)
    update_drainage_info(model, divides_x)
    return topo_gridx, partial_melt_gridy, divides_x
end

""" Sorted-x equivalent of `calculate_melt_drainage_divides!`.

Identical inputs, outputs, and side effects as the original; only the call to
`calculate_top_of_mantle_partial_melt_domain` is swapped for its `_sorted`
variant, which uses sorted-x binary-search bracketing instead of a linear
x-range scan over filtered markers.
"""
function calculate_melt_drainage_divides_sorted!(
    model::ModelData
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    topo_gridx = model.topography.arrays.gridt.array[1, :]
    partial_melt_gridy = calculate_top_of_mantle_partial_melt_domain_sorted(
        model, topo_gridx)
    divides_x = calculate_drainage_divides(topo_gridx, partial_melt_gridy)
    update_drainage_info(model, divides_x)
    return topo_gridx, partial_melt_gridy, divides_x
end

""" Update drainage information.

The drainage divide arrays stored in the main data object are updated with the updated 
drainage divides.

# Updated Arrays
- model.melting.arrays.extraction
    - xstart_drainage : Vector{Float64}
        X-start locations of drainage divides (meters)
    - xend_drainage : Vector{Float64}
        X-end locations of drainage divides (meters)
"""
function update_drainage_info(
    model::ModelData,
    divides_x::Vector{Float64}
)::Nothing
    copy_old_number_of_drainage_basins(model)
    ndrainage_basin, ndivides = update_current_number_of_drainage_basins(model, divides_x)
    xstart_drainage   = model.melting.arrays.extraction.xstart_drainage.array
    xend_drainage     = model.melting.arrays.extraction.xend_drainage.array
    xstart_drainage_o = model.melting.arrays.extraction.xstart_drainage_o.array
    xend_drainage_o   = model.melting.arrays.extraction.xend_drainage_o.array
    nmax = length(xstart_drainage)
    if ndivides > nmax
        println(">> Maximum number of drainage divides: ", nmax)
        println(
            "!!! WARNING !!! Number of drainage divides exceeds maximum " *
            "number of drainage divides. The excess will be ignored."
        )
        ndivides = nmax
    end
    # Save the old drainage divide arrays
    xstart_drainage_o[1:nmax] .= xstart_drainage[1:nmax]
    xend_drainage_o[1:nmax] .= xend_drainage[1:nmax]
    # Clear the current drainage divide arrays
    xstart_drainage[1:nmax] .= 0.0
    xend_drainage[1:nmax] .= 0.0
    # Update the current drainage divide arrays
    xstart_drainage[1:ndrainage_basin] .= divides_x[1:ndrainage_basin]
    xend_drainage[1:ndrainage_basin] .= divides_x[2:ndrainage_basin+1]
    #for i in 1:ndrainage_basin
    #    xend_drainage[i] = divides_x[i+1]
    #end
    print_info = false
    if print_info
        print_drainage_info(ndivides, ndrainage_basin, xstart_drainage, xend_drainage)
    end
    return nothing
end

function copy_old_number_of_drainage_basins(model::ModelData)::Nothing
    ndrainage_basin = model.melting.parameters.extraction.ndrainage_basin.value
    model.melting.parameters.extraction.ndrainage_basin_old.value = ndrainage_basin
    return nothing
end

function update_current_number_of_drainage_basins(
    model::ModelData, 
    divides_x::Vector{Float64}
)::Tuple{Int64, Int64}
    ndivides = length(divides_x)
    ndrainage_basin = ndivides - 1
    model.melting.parameters.extraction.ndrainage_basin.value = ndrainage_basin
    return ndrainage_basin, ndivides
end

function print_drainage_info(
    ndivides::Int64,
    ndrainage_basin::Int64,
    xstart_drainage::Vector{Float64},
    xend_drainage::Vector{Float64}
)::Nothing
    print_info("Number of drainage divides: $(ndivides)", level=1)
    print_info("Number of drainage basins: $(ndrainage_basin)", level=1)
    for i in 1:ndrainage_basin
        print_info(
            "Drainage basin $(i) xstart $(xstart_drainage[i]) xend $(xend_drainage[i])", level=2
        )
    end
    print_info("", level=1)
    return nothing
end

function calculate_top_of_mantle_partial_melt_domain(
    model::ModelData,
    topo_gridx::Vector{Float64}
)::Vector{Float64}
    smoothing_radius = model.melting.parameters.extraction.smoothing_radius_drainage.value
    ysize = model.grids.parameters.geometry.ysize.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    matid_types = model.materials.dicts.matid_types
    matids_mantle_partial_melt = matid_types["UltramaficMantlePartiallyMolten"]

    dx_grid = topo_gridx[2] - topo_gridx[1]

    mxstep = model.markers.parameters.distribution.mxstep.value
    marker_search_factor = model.topography.parameters.topo_grid.marker_search_factor.value
    search_radius = calculate_search_radius(mxstep, topo_gridx, marker_search_factor)

    # Do not use smoothing for this case since smoothing will be done after
    # setting zero values to model base
    (
        top_mantle_partial_melt, _bottom_mantle_partial_melt
    ) = calculate_top_and_bottom_of_layer_opt(
        matids_mantle_partial_melt, marker_matid, marker_x, marker_y,
        topo_gridx, search_radius; use_smoothing=false
    )
    set_zero_values_to_model_base(top_mantle_partial_melt, ysize)

    nsmooth = floor(Int, smoothing_radius / dx_grid)

    debug = false
    if debug
        println(">> nsmooth for partial melt domain: ", nsmooth)
        println(
            ">> Min and max top_mantle_partial_melt (before smoothing): ",
            minimum(top_mantle_partial_melt), " ", maximum(top_mantle_partial_melt)
        )
    end

    top_mantle_partial_melt = smooth_surface(top_mantle_partial_melt; nsmooth=nsmooth)

    if debug
        println(
            ">> Min and max top_mantle_partial_melt (after smoothing): ",
            minimum(top_mantle_partial_melt), " ", maximum(top_mantle_partial_melt)
        )
    end

    return top_mantle_partial_melt
end

""" Sorted-x equivalent of `calculate_top_of_mantle_partial_melt_domain`.

Identical inputs and outputs as the original; the only difference is that the
underlying layer scan uses `calculate_top_and_bottom_of_layer_opt_sorted`,
which brackets the marker x-range per grid column via binary search instead
of linearly scanning every filtered marker.
"""
function calculate_top_of_mantle_partial_melt_domain_sorted(
    model::ModelData,
    topo_gridx::Vector{Float64}
)::Vector{Float64}
    smoothing_radius = model.melting.parameters.extraction.smoothing_radius_drainage.value
    ysize = model.grids.parameters.geometry.ysize.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    matid_types = model.materials.dicts.matid_types
    matids_mantle_partial_melt = matid_types["UltramaficMantlePartiallyMolten"]

    dx_grid = topo_gridx[2] - topo_gridx[1]

    mxstep = model.markers.parameters.distribution.mxstep.value
    marker_search_factor = model.topography.parameters.topo_grid.marker_search_factor.value
    search_radius = calculate_search_radius(mxstep, topo_gridx, marker_search_factor)

    (
        top_mantle_partial_melt, _bottom_mantle_partial_melt
    ) = calculate_top_and_bottom_of_layer_opt_sorted(
        matids_mantle_partial_melt, marker_matid, marker_x, marker_y,
        topo_gridx, search_radius; use_smoothing=false
    )
    set_zero_values_to_model_base(top_mantle_partial_melt, ysize)

    nsmooth = floor(Int, smoothing_radius / dx_grid)
    top_mantle_partial_melt = smooth_surface(top_mantle_partial_melt; nsmooth=nsmooth)

    return top_mantle_partial_melt
end

function set_zero_values_to_model_base(array::Vector{Float64}, ysize::Float64)::Nothing
    n = length(array)
    Threads.@threads for i in 1:n
        if array[i] == 0.0
            array[i] = ysize
        end
    end
    return nothing
end

""" Calculate drainage divides.

This function calculates drainage divides for a curve defined by the input x and y 
coordinates. The drainage divides separate the curve into regions where the curve 
slopes downward in opposite directions.

The y-coordinate is positive downward. Therefore, the slope is positive when the 
y-coordinate increases with increasing x-coordinate.

Parameters
----------
topo_gridx : Vector{Float64}
    Grid x-coordinates (meters)
topo_gridy : Vector{Float64}
    Grid y-coordinates (meters). Note that y increases with depth

Returns
-------
divides_x : Vector{Float64}
    X-locations of drainage divides (meters)
"""
function calculate_drainage_divides(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64}
)::Vector{Float64}
    divides_x_tmp = Vector{Float64}(undef, length(topo_gridx))
    # Add left end point to divides_x_tmp
    divides_x_tmp[1] = topo_gridx[1]
    icount = 2
    # Add internal divides to divides_x_tmp
    for i in 2:length(topo_gridx)-1
        left_slope = (topo_gridy[i] - topo_gridy[i-1]) / (topo_gridx[i] - topo_gridx[i-1])
        right_slope = (topo_gridy[i+1] - topo_gridy[i]) / (topo_gridx[i+1] - topo_gridx[i])
        
        if left_slope > 0.0 && right_slope < 0.0
            divides_x_tmp[icount] = topo_gridx[i]
            icount += 1
        elseif left_slope > 0.0 && right_slope == 0.0
            x_start = topo_gridx[i]
            right_slope_b = 0.0  # Initialize before the loop
            x_end = 0.0 # Initialize before the loop
            for m in i+1:length(topo_gridx)-1
                right_slope_b = (topo_gridy[m+1] - topo_gridy[m]) / (topo_gridx[m+1] - topo_gridx[m])
                if right_slope_b != 0.0
                    x_end = topo_gridx[m]
                    break
                end
            end
            if right_slope_b < 0.0
                divides_x_tmp[icount] = (x_start + x_end) / 2
                icount += 1
            end
        end
    end
    # Add right end point to divides_x_tmp
    divides_x_tmp[icount] = topo_gridx[end]
    icount += 1
    # Clean up divides_x_tmp
    divides_x = Vector{Float64}(undef, icount-1)
    for i in 1:icount-1
        divides_x[i] = divides_x_tmp[i]
    end
    return divides_x
end

""" Redistribute accumulated extrusion volumes from an old set of drainage
basins to a new set, conserving total volume.

Each old basin's volume is assigned to exactly ONE new basin: the new basin
that contains the old basin's midpoint. New basins tile the model domain
contiguously and share exact boundaries, so a midpoint that lands on a shared
boundary is resolved with a half-open rule `[xstart, xend)` — i.e. it is
assigned to the basin on the right. (The previous inclusive-both-ends test
`xstart <= xmid_o <= xend` matched both neighbouring basins for such a midpoint
and added the volume twice, doubling it.) If a midpoint falls outside every new
basin — e.g. the new basins do not tile the domain — its volume is assigned to
the nearest basin so the total is still conserved.

Operates on the first `ndrainage_basin` / `ndrainage_basin_o` entries of the
arrays. `extrusion_volumes` is fully zeroed and then filled; it must be a
distinct array from `extrusion_volumes_o` (the old volumes).
"""
function redistribute_extrusion_volumes!(
    extrusion_volumes::Vector{Float64},
    extrusion_volumes_o::Vector{Float64},
    xstart_drainage::Vector{Float64},
    xend_drainage::Vector{Float64},
    xstart_drainage_o::Vector{Float64},
    xend_drainage_o::Vector{Float64},
    ndrainage_basin::Int,
    ndrainage_basin_o::Int
)::Nothing
    fill!(extrusion_volumes, 0.0)
    for j in 1:ndrainage_basin_o
        xmid_o = (xstart_drainage_o[j] + xend_drainage_o[j]) / 2
        inew = 0
        for i in 1:ndrainage_basin
            xstart = xstart_drainage[i]
            xend = xend_drainage[i]
            # Half-open [xstart, xend) assigns a boundary midpoint to the basin
            # on the right; the last basin also owns the domain's right edge.
            if (xstart <= xmid_o < xend) || (i == ndrainage_basin && xmid_o <= xend)
                inew = i
                break
            end
        end
        if inew == 0
            # Midpoint outside every new basin: assign to the nearest basin so
            # no volume is lost.
            best_distance = Inf
            for i in 1:ndrainage_basin
                xmid = (xstart_drainage[i] + xend_drainage[i]) / 2
                distance = abs(xmid_o - xmid)
                if distance < best_distance
                    best_distance = distance
                    inew = i
                end
            end
        end
        extrusion_volumes[inew] += extrusion_volumes_o[j]
    end
    return nothing
end

function redistribute_extrusion_volumes_to_new_drainage_basins!(model::ModelData)::Nothing
    ndrainage_basin_o = model.melting.parameters.extraction.ndrainage_basin_old.value
    ndrainage_basin = model.melting.parameters.extraction.ndrainage_basin.value

    xstart_drainage = model.melting.arrays.extraction.xstart_drainage.array
    xend_drainage = model.melting.arrays.extraction.xend_drainage.array
    xstart_drainage_o = model.melting.arrays.extraction.xstart_drainage_o.array
    xend_drainage_o = model.melting.arrays.extraction.xend_drainage_o.array

    extrusion_volumes = model.melting.arrays.extraction.extrusion_volumes.array
    extrusion_volumes_o = copy(extrusion_volumes)

    # Calculate current total extrusion volume from all drainage basins
    total_extrusion_volume_current = sum(extrusion_volumes)

    redistribute_extrusion_volumes!(
        extrusion_volumes, extrusion_volumes_o,
        xstart_drainage, xend_drainage,
        xstart_drainage_o, xend_drainage_o,
        ndrainage_basin, ndrainage_basin_o
    )

    # Calculate new total extrusion volume from all drainage basins
    total_extrusion_volume_redistributed = sum(extrusion_volumes)

    # Conservation is guaranteed by construction; this is a loud, non-fatal
    # guard against future regressions (uses @warn so it is never suppressed by
    # PRINT_SETTINGS, unlike the EarthBox print helpers).
    if !isapprox(total_extrusion_volume_current, total_extrusion_volume_redistributed)
        @warn(
            "Total extrusion volume is not conserved after redistributing " *
            "extrusion volumes to new drainage basins.",
            total_extrusion_volume_current,
            total_extrusion_volume_redistributed,
            difference = total_extrusion_volume_current - total_extrusion_volume_redistributed,
        )
    else
        print_info(
            "Total extrusion volume before and after redistribution: " *
            "$(total_extrusion_volume_current) $(total_extrusion_volume_redistributed)",
            level=2
        )
    end

    return nothing
end

end # module 