module CompactionCorrection

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelStructureManager.TopAndBottom: calculate_layer_thickness
import EarthBox.ModelStructureManager.TopAndBottom: calculate_search_radius
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm_opt
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm
import EarthBox.DataStructures: SedimentTransportParameters
import ..CompactionTools: compact_or_decompact_layer
import ..MarkerCompaction: compact_sediment_and_advect_markers

""" Apply compaction correction to transport topography and markers.

# Arguments
- `model`: Model data
- `sediment_transport_parameters`: Sediment transport parameters
- `topo_gridx`: Topography grid x-coordinates (meters)
- `topo_gridy_initial`: Initial topography grid y-coordinates (meters)
- `topo_gridy_new`: New topography grid y-coordinates (meters)
- `sediment_thickness_initial`: Initial sediment thickness (meters)
- `markers_indices_sedimentary_basin`: Indices of markers in sedimentary basin
- `markers_indices_sticky`: Indices of sticky markers
- `x_sorted_marker_indices_sedimentary_basin`: X-sorted indices of markers in 
    sedimentary basin
- `x_sorted_marker_indices_sticky`: X-sorted indices of sticky markers

# Returns
- `topo_gridy_corrected`: Corrected topography grid y-coordinates (meters)
- `total_sed_thickness_corrected`: Corrected total sediment thickness (meters)
- `new_thickness_decompacted`: Decompacted new sediment thickness (meters)
"""
function apply_compaction_correction_for_topography_and_markers(
    model::ModelData,
    sediment_transport_parameters::SedimentTransportParameters,
    topo_gridx::Vector{Float64},
    topo_gridy_initial::Vector{Float64},
    topo_gridy_new::Vector{Float64},
    sediment_thickness_initial::Vector{Float64},
    markers_indices_sedimentary_basin::Vector{Int64},
    markers_indices_sticky::Vector{Int64},
    x_sorted_marker_indices_sedimentary_basin::Vector{Int64},
    x_sorted_marker_indices_sticky::Vector{Int64}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    markers_x = model.markers.arrays.location.marker_x.array
    markers_y = model.markers.arrays.location.marker_y.array
    marker_porosity_initial = model.markers.arrays.material.marker_porosity_initial.array
    marker_decay_depth = model.markers.arrays.material.marker_decay_depth.array
    markers_max_burial_depth = model.markers.arrays.material.marker_max_burial_depth.array
    marker_search_factor = model.topography.parameters.topo_grid.marker_search_factor.value
    nsmooth_top_bottom = model.topography.parameters.topo_grid.nsmooth_top_bottom.value
    mxstep = model.markers.parameters.distribution.mxstep.value
   
    search_radius = calculate_search_radius(mxstep, topo_gridx, marker_search_factor)
    top_sticky, bottom_sticky = calculate_top_and_bottom_of_swarm_opt(
        x_sorted_marker_indices_sticky,
        markers_x, markers_y, topo_gridx,
        search_radius,
        use_smoothing=true,
        nsmooth=nsmooth_top_bottom
    )
    
    sticky_thickness_initial = bottom_sticky .- top_sticky
    
    porosity_initial_transport = sediment_transport_parameters.porosity_initial
    decay_depth_transport = 1.0/sediment_transport_parameters.depth_decay_term
    
    (
        topo_gridy_corrected, total_sed_thickness_corrected, new_thickness_decompacted
    ) = decompact_new_sediment_and_compact_markers(
            porosity_initial_transport,
            decay_depth_transport,
            topo_gridx,
            topo_gridy_initial,
            topo_gridy_new,
            markers_x,
            markers_y,
            marker_porosity_initial,
            marker_decay_depth,
            markers_max_burial_depth,
            sediment_thickness_initial,
            sticky_thickness_initial,
            markers_indices_sedimentary_basin,
            markers_indices_sticky,
            x_sorted_marker_indices_sedimentary_basin,
            search_radius,
            nsmooth_top_bottom
        )
    
    return topo_gridy_corrected, total_sed_thickness_corrected, new_thickness_decompacted
end

"""
Apply compaction correction to transport topography and markers.

# Updated Arrays
- `markers_y`: Array of marker y-coordinates (meters)
- `markers_max_burial_depth`: Array of marker maximum burial depths (meters)

# Algorithm
1. Calculate new sediment thickness with full compaction by taking the
   difference between the initial and new topography grids.
2. Decompact the new sediment thickness.
3. Compact the initial sediment thickness using the decompacted new sediment
   thickness and apply compaction to sediment markers and sticky markers.
4. Calculate the total sediment thickness before and after compaction.
5. Correct the new topography grid by subtracting the difference between
   the total sediment thickness before and after compaction:

    topo_gridy_corrected = (
        topo_gridy_new
        - (total_sed_thickness_corrected - total_sed_thickness_uncorrected)
        )

# Arguments
- `porosity_initial_transport::Float64`: Initial porosity of newly deposited sediment (fraction)
- `decay_depth_transport::Float64`: Depth decay term of deposited sediment (meters)
- `topo_gridx::Vector{Float64}`: X coordinates of topography grid (meters)
- `topo_gridy_initial::Vector{Float64}`: Initial topography grid y-coordinates (meters)
- `topo_gridy_from_transport::Vector{Float64}`: New topography grid y-coordinates (meters)
- `markers_x::Vector{Float64}`: Marker x-coordinates (meters)
- `markers_y::Vector{Float64}`: Marker y-coordinates (meters)
- `marker_porosity_initial::Vector{Float32}`: Initial marker porosity (fraction)
- `marker_decay_depth::Vector{Float32}`: Marker porosity decay depth (meters)
- `markers_max_burial_depth::Vector{Float32}`: Maximum burial depth of markers (meters)
- `sediment_thickness_initial::Vector{Float64}`: Initial sediment thickness (meters)
- `sticky_thickness_initial::Vector{Float64}`: Initial sticky thickness (meters)
- `markers_indices_sedimentary_basin::Vector{Int64}`: Indices of sedimentary basin markers
- `markers_indices_sticky::Vector{Int64}`: Indices of sticky markers
- `x_sorted_marker_indices_sedimentary_basin::Vector{Int64}`: X-sorted sedimentary basin marker indices
- `search_radius::Float64`: Search radius for markers (meters)
- `nsmooth_top_bottom::Int64=2`: Number of smoothing iterations

# Returns
- `topo_gridy_corrected::Vector{Float64}`: Corrected topography y-coordinates (meters)
- `total_sed_thickness_corrected::Vector{Float64}`: Total corrected sediment thickness (meters)
- `new_thickness_decompacted::Vector{Float64}`: Decompacted new sediment thickness (meters)
"""
function decompact_new_sediment_and_compact_markers(
    porosity_initial_transport::Float64,
    decay_depth_transport::Float64,
    topo_gridx::Vector{Float64},
    topo_gridy_initial::Vector{Float64},
    topo_gridy_from_transport::Vector{Float64},
    markers_x::Vector{Float64},
    markers_y::Vector{Float64},
    marker_porosity_initial::Vector{Float32},
    marker_decay_depth::Vector{Float32},
    markers_max_burial_depth::Vector{Float32},
    sediment_thickness_initial::Vector{Float64},
    sticky_thickness_initial::Vector{Float64},
    markers_indices_sedimentary_basin::Vector{Int64},
    markers_indices_sticky::Vector{Int64},
    x_sorted_marker_indices_sedimentary_basin::Vector{Int64},
    search_radius::Float64,
    nsmooth_top_bottom::Int64=2
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    new_thickness_compacted = calculate_new_sediment_thickness_full_compaction(
        topo_gridy_initial, topo_gridy_from_transport
    )
    new_thickness_decompacted = decompact_transport_thickness_opt(
        porosity_initial_transport, decay_depth_transport,
        new_thickness_compacted
    )
    sediment_thickness_initial_compacted = compact_sediment_and_advect_markers(
        topo_gridx,
        topo_gridy_initial,
        sediment_thickness_initial,
        sticky_thickness_initial,
        new_thickness_decompacted,
        marker_porosity_initial,
        marker_decay_depth,
        markers_max_burial_depth,
        markers_indices_sedimentary_basin,
        markers_indices_sticky,
        x_sorted_marker_indices_sedimentary_basin,
        markers_x,
        markers_y,
        search_radius,
        nsmooth_top_bottom;
        default_porosity_initial=porosity_initial_transport,
        default_decay_depth=decay_depth_transport
    )
    total_sed_thickness_uncorrected = sediment_thickness_initial .+ new_thickness_compacted
    total_sed_thickness_corrected = sediment_thickness_initial_compacted .+ new_thickness_decompacted
    topo_gridy_corrected = topo_gridy_from_transport .- 
        (total_sed_thickness_corrected .- total_sed_thickness_uncorrected)
    return (
        topo_gridy_corrected,
        total_sed_thickness_corrected,
        new_thickness_decompacted
    )
end

""" Apply compaction correction.

# Arguments
- `sediment_transport_parameters`: Sediment transport parameters
- `topo_gridy_initial`: Initial topography grid y-coordinates (meters)
- `topo_gridy_new`: New topography grid y-coordinates (meters)
- `sediment_thickness_initial`: Initial sediment thickness (meters)

# Returns
- `topo_gridy_corrected`: Corrected topography grid y-coordinates (meters)
- `total_sed_thickness_corrected`: Corrected total sediment thickness (meters)
- `new_thickness_decompacted`: Decompacted new sediment thickness (meters)
"""
function apply_compaction_correction(
    sediment_transport_parameters::SedimentTransportParameters,
    topo_gridy_initial::Vector{Float64},
    topo_gridy_new::Vector{Float64},
    sediment_thickness_initial::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    new_thickness_compacted = calculate_new_sediment_thickness_full_compaction(
        topo_gridy_initial, topo_gridy_new
    )
    
    new_thickness_decompacted = decompact_transport_thickness(
        sediment_transport_parameters, new_thickness_compacted
    )
    
    sediment_thickness_initial_compacted = compact_sediment_thickness_initial(
        sediment_transport_parameters,
        sediment_thickness_initial,
        new_thickness_decompacted
    )
    
    total_sed_thickness_uncorrected = sediment_thickness_initial .+ new_thickness_compacted
    total_sed_thickness_corrected = sediment_thickness_initial_compacted .+ 
        new_thickness_decompacted
    
    topo_gridy_corrected = topo_gridy_new .- 
        (total_sed_thickness_corrected .- total_sed_thickness_uncorrected)
    
    return topo_gridy_corrected, total_sed_thickness_corrected, new_thickness_decompacted
end

""" Calculate new sediment thickness with decompaction.

# Arguments
- `sediment_transport_parameters`: Sediment transport parameters
- `sediment_thickness_transport_compacted`: Compacted sediment thickness (meters)

# Returns
- `sediment_thickness_transport_decompacted`: Decompacted sediment thickness (meters)
"""
function decompact_transport_thickness(
    sediment_transport_parameters::SedimentTransportParameters,
    sediment_thickness_transport_compacted::Vector{Float64}
)::Vector{Float64}
    porosity_initial = sediment_transport_parameters.porosity_initial
    depth_decay_term = sediment_transport_parameters.depth_decay_term
    xnum = length(sediment_thickness_transport_compacted)
    # Use a very deep depth since the sediment thickness is fully compacted
    top_initial_at_15km = fill(15_000.0, xnum)
    # Decompact to a zero depth below mudline
    top_new_at_zero = zeros(xnum)
    sediment_thickness_transport_decompacted = compact_or_decompact_layer(
        porosity_initial,
        depth_decay_term,
        top_initial_at_15km,
        sediment_thickness_transport_compacted,
        top_new_at_zero
    )
    return sediment_thickness_transport_decompacted
end

""" Calculate new sediment thickness with decompaction.

# Arguments
- `porosity_initial`: Initial porosity (fraction)
- `depth_decay`: Depth decay depth (meters)
- `sediment_thickness_transport_compacted`: Compacted sediment thickness (meters)

# Returns
- `sediment_thickness_transport_decompacted`: Decompacted sediment thickness (meters)
"""
function decompact_transport_thickness_opt(
    porosity_initial::Float64,
    depth_decay::Float64,
    sediment_thickness_transport_compacted::Vector{Float64}
)::Vector{Float64}
    xnum = length(sediment_thickness_transport_compacted)
    # Use a very deep depth since the sediment thickness is fully compacted
    top_initial_at_15km = 15_000.0
    # Decompact to a zero depth below mudline
    top_new_at_zero = 0.0
    depth_decay_term = 1.0/depth_decay
    sediment_thickness_transport_decompacted = compact_or_decompact_layer(
        porosity_initial,
        depth_decay_term,
        top_initial_at_15km,
        sediment_thickness_transport_compacted,
        top_new_at_zero
    )
    
    return sediment_thickness_transport_decompacted
end

""" Calculate new sediment thickness with full compaction.

# Arguments
- `topo_gridy_initial`: Initial topography grid y-coordinates (meters)
- `topo_gridy_new`: New topography grid y-coordinates (meters)

# Returns
- `sediment_thickness_compacted_new`: Compacted new sediment thickness (meters)
"""
function calculate_new_sediment_thickness_full_compaction(
    topo_gridy_initial::Vector{Float64},
    topo_gridy_new::Vector{Float64}
)::Vector{Float64}
    topo_delta_y = calculate_change_in_topography(topo_gridy_initial, topo_gridy_new)
    sediment_thickness_compacted_new = max.(0.0, topo_delta_y)
    return sediment_thickness_compacted_new
end

""" Calculate change in topography.

# Arguments
- `topo_gridy_initial`: Initial topography grid y-coordinates (meters)
- `topo_gridy_new`: New topography grid y-coordinates (meters)

# Returns
- `topo_delta_y`: Change in topography (meters). A positive delta y indicates 
    deposition and a negative delta y indicates erosion.
"""
function calculate_change_in_topography(
    topo_gridy_initial::Vector{Float64},
    topo_gridy_new::Vector{Float64}
)::Vector{Float64}
    topo_delta_y = topo_gridy_initial - topo_gridy_new
    return topo_delta_y
end

""" Decompact sediment thickness from markers.

# Arguments
- `sediment_transport_parameters`: Sediment transport parameters
- `sediment_thickness_initial`: Initial sediment thickness (meters)
- `sediment_thickness_transport_decompacted`: Decompacted transport sediment 
    thickness (meters)

# Returns
- `sediment_thickness_original_compacted`: Compacted original sediment thickness 
    (meters)
"""
function compact_sediment_thickness_initial(
    sediment_transport_parameters::SedimentTransportParameters,
    sediment_thickness_initial::Vector{Float64},
    sediment_thickness_transport_decompacted::Vector{Float64}
)::Vector{Float64}
    porosity_initial = sediment_transport_parameters.porosity_initial
    depth_decay_term = sediment_transport_parameters.depth_decay_term
    
    xnum = length(sediment_thickness_transport_decompacted)
    top_initial = zeros(xnum)
    
    sediment_thickness_original_compacted = compact_or_decompact_layer(
        porosity_initial,
        depth_decay_term,
        top_initial,  # Initial submud top
        sediment_thickness_initial,  # Initial thickness
        sediment_thickness_transport_decompacted  # New submud top
    )
    
    return sediment_thickness_original_compacted
end

end # module 