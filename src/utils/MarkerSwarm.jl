module MarkerSwarm

""" Calculate dimensions of a swarm of markers with a given set of material ids.

# Arguments
- marker_x::Vector{Float64}
    - X-coordinates of markers.
- marker_y::Vector{Float64}
    - Y-coordinates of markers.
- marker_matid::Vector{Int16}
    - Material IDs of markers.
- target_matids::Vector{Int16}
    - Target material IDs to search for.
- xstart::Float64
    - Starting x-coordinate of drainage basin.
- xend::Float64
    - Ending x-coordinate of drainage basin.

# Returns
- xmid_target::Float64
    - X-location of the middle of target zone.
- ytop_target::Float64
    - Minimum y-location of the target zone.
- width_target::Float64
    - Width of target zone.
"""
function calculate_dimensions_general(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    target_matids::Vector{Int16},
    xstart::Float64 = -1e39,
    xend::Float64 = 1e39
)::Tuple{Float64, Float64, Float64}
    # Find dimensions of target zone
    xmin_mol = 1e32
    xmax_mol = -1e32
    ymin_mol = 1e32
    ymax_mol = -1e32
    marknum = length(marker_x)
    
    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        if xstart <= x_marker <= xend
            y_marker = marker_y[imarker]
            matid = marker_matid[imarker]
            if matid in target_matids
                if y_marker < ymin_mol
                    ymin_mol = y_marker
                end
                if y_marker > ymax_mol
                    ymax_mol = y_marker
                end
                if x_marker < xmin_mol
                    xmin_mol = x_marker
                end
                if x_marker > xmax_mol
                    xmax_mol = x_marker
                end
            end
        end
    end
    
    xmid_target = (xmax_mol + xmin_mol)/2.0
    ytop_target = ymin_mol
    width_target = xmax_mol - xmin_mol
    return xmid_target, ytop_target, width_target
end

end # module