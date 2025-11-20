module MagmaBodyTest

import EarthBox.MeltModel.Extraction: PartiallyMoltenZone
import EarthBox.FindShallowest: find_shallowest_partially_molten_mantle_marker_opt
import Plots

function run_test()
    (
        marker_x, marker_y, 
        marker_matid, partial_melt_flags, mantle_melting_mat_ids, 
        nmarkers_partial_melt, nmarkers_magma, matid_magma
    ) = make_inputs()

    (
        layer_counts, marker_indices
    ) = PartiallyMoltenZone.construct_layered_partially_molten_arrays(
        marker_x, marker_y, marker_matid, nmarkers_partial_melt,
        mantle_melting_mat_ids, partial_melt_flags;
        nlayers = 10
    )

    for _ in 1:nmarkers_magma
        (
            imarker_shallow,
            _yshallow
        ) = find_shallowest_partially_molten_mantle_marker_opt(
            marker_matid,
            marker_y,
            mantle_melting_mat_ids,
            layer_counts,
            marker_indices
        )

        if imarker_shallow != -999
            marker_matid[imarker_shallow] = matid_magma
        end
    end

    plot_results(marker_x, marker_y, marker_matid)
end

function plot_results(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16}
)::Nothing
    dpi = 150
    figsize = (10, 10)
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)
    p = Plots.scatter(
        marker_x ./ 1000.0, 
        marker_y ./ 1000.0, 
        marker_z=marker_matid,
        clims=(1.0, 3.0),
        color=:viridis,
        colorbar_title="marker_matid",
        xlabel="X Axis (km)",
        ylabel="Y Axis (km)",
        title="Magma body use_optimized",
        yflip=true,
        aspect_ratio=:auto,
        size=figsize_pixels,
        markershape=:rect,
        markersize=4.0,
        markerstrokewidth=0.0,
        markerstrokecolor=:transparent,
        )
    Plots.savefig(p, "magma_body.png")
    return nothing
end

""" Make inputs.

# Arguments
- `nmarkers_x::Int`: Number of markers in x direction
- `nmarkers_y::Int`: Number of markers in y direction
- `dx::Float64`: Marker spacing in x direction
- `dy::Float64`: Marker spacing in y direction
- `matid_background::Int`: Background material id
- `matid_pm::Int`: Partial melt material id
- `matid_magma::Int`: Magma material id
- `magma_fraction::Float64`: Fraction of partially molten markers that will be 
                           converted into magma

# Returns
- `marker_x::Vector{Float64}`: Marker x coordinates
- `marker_y::Vector{Float64}`: Marker y coordinates
- `marker_matid::Vector{Int16}`: Marker material ids
- `partial_melt_flags::Vector{Int64}`: Array of marker indices that are partially molten
- `mantle_melting_mat_ids::Vector{Int64}`: Material ids of mantle melting markers
- `nmarkers_partial_melt::Int`: Number of markers that are partially molten
- `nmarkers_magma::Int`: Number of markers that are magma
"""
function make_inputs(;
    nmarkers_x::Int=101,
    nmarkers_y::Int=101,
    dx::Float64=1.0,
    dy::Float64=1.0,
    matid_background::Int16=Int16(1),
    matid_pm::Int16=Int16(2),
    matid_magma::Int16=Int16(3),
    magma_fraction::Float64=0.3
)::Tuple{
    Vector{Float64}, Vector{Float64},
    Vector{Int16}, Vector{Int64}, Vector{Int16}, 
    Int, Int, Int16
}
    (
        marker_x, marker_y, marker_matid
    ) = make_markers(
        nmarkers_x, nmarkers_y, matid_background, dx=dx, dy=dy)

    partial_melt_flags = make_partial_melt_flags(
        marker_x, marker_y, marker_matid, matid_pm)

    nmarkers_partial_melt = length(partial_melt_flags)
    nmarkers_magma = floor(Int, magma_fraction * nmarkers_partial_melt)
    mantle_melting_mat_ids = [Int16(2)]

    return (
        marker_x, marker_y, 
        marker_matid, partial_melt_flags, mantle_melting_mat_ids, 
        nmarkers_partial_melt, nmarkers_magma, matid_magma
    )
end

""" Make markers and initialized material id's using background material.

# Arguments
- `nmarkers_x::Int`: Number of markers in x direction
- `nmarkers_y::Int`: Number of markers in y direction
- `matid_background::Int16`: Background material id
- `dx::Float64`: Marker spacing in x direction
- `dy::Float64`: Marker spacing in y direction

# Returns
- `marker_x::Vector{Float64}`: Marker x coordinates
- `marker_y::Vector{Float64}`: Marker y coordinates
- `marker_matid::Vector{Int16}`: Marker material ids
"""
function make_markers(
    nmarkers_x::Int,
    nmarkers_y::Int,
    matid_background::Int16;
    dx::Float64=1.0,
    dy::Float64=1.0
)::Tuple{
    Vector{Float64},
    Vector{Float64},
    Vector{Int16}
}
    nmarkers = nmarkers_x * nmarkers_y
    marker_x = zeros(Float64, nmarkers)
    marker_y = zeros(Float64, nmarkers)
    marker_matid = zeros(Int16, nmarkers)
    icount = 1
    for j in 1:nmarkers_y
        for i in 1:nmarkers_x
            x_marker = Float64(i-1) * dx
            y_marker = Float64(j-1) * dy
            matid = matid_background
            marker_x[icount] = x_marker
            marker_y[icount] = y_marker
            marker_matid[icount] = matid
            icount += 1
        end
    end
    return marker_x, marker_y, marker_matid
end

""" Make partial melt flags.

# Arguments
- `marker_x::Vector{Float64}`: Marker x coordinates
- `marker_y::Vector{Float64}`: Marker y coordinates
- `marker_matid::Vector{Int16}`: Marker material ids
- `matid_pm::Int16`: Partial melt material id
- `left_edge_molten_box_fraction::Float64`: Location of left edge of molten box in 
                                          fraction of max x
- `right_edge_molten_box_fraction::Float64`: Location of right edge of molten box in 
                                           fraction of max x
- `top_edge_molten_box_fraction::Float64`: Location of top edge of molten box in 
                                         fraction of max y
- `bottom_edge_molten_box_fraction::Float64`: Location of bottom edge of molten box in 
                                            fraction of max y

# Returns
- `partial_melt_flags::Vector{Int64}`: Array of marker indices that are partially molten
"""
function make_partial_melt_flags(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    matid_pm::Int16;
    left_edge_molten_box_fraction::Float64=0.25,
    right_edge_molten_box_fraction::Float64=0.75,
    top_edge_molten_box_fraction::Float64=0.1,
    bottom_edge_molten_box_fraction::Float64=0.35
)::Vector{Int64}
    xpm_o = maximum(marker_x) * left_edge_molten_box_fraction
    xpm_f = maximum(marker_x) * right_edge_molten_box_fraction
    ypm_o = maximum(marker_y) * top_edge_molten_box_fraction
    ypm_f = maximum(marker_y) * bottom_edge_molten_box_fraction

    println(">> xpm_o: ", xpm_o)
    println(">> xpm_f: ", xpm_f)
    println(">> ypm_o: ", ypm_o)
    println(">> ypm_f: ", ypm_f)

    nmarkers = length(marker_x)
    partial_melt_flags_tmp = zeros(Int64, nmarkers)
    icount = 1
    for imarker in 1:nmarkers
        if marker_x[imarker] >= xpm_o && marker_x[imarker] <= xpm_f
            if marker_y[imarker] >= ypm_o && marker_y[imarker] <= ypm_f
                marker_matid[imarker] = matid_pm
                partial_melt_flags_tmp[icount] = imarker
                icount += 1
            end
        end
    end
    partial_melt_flags = zeros(Int64, icount-1)
    for i in 1:(icount-1)
        partial_melt_flags[i] = partial_melt_flags_tmp[i]
    end
    println(">> nmarkers_partial_melt: ", icount-1)
    #println(">> partial_melt_flags: ", partial_melt_flags)
    return partial_melt_flags
end

end # module