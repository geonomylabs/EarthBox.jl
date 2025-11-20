module ApplyCompaction

import EarthBox.ModelDataContainer: ModelData
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.EBCopy: copy_array_1d!
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm_opt
import ..MarkerCompaction: calculate_swarm_indices_for_sediment_and_sticky
import ..MarkerCompaction: calculate_x_sorted_swarm_indices
import ..CompactionCorrection: apply_compaction_correction_for_topography_and_markers

""" Apply the compaction model.

# Arguments
- model: Union{ModelData, Nothing}
    The model data object.

- topo_gridx: Vector{Float64}
    The x-coordinates of the topography grid (meters).

- topo_gridy: Vector{Float64}
    The y-coordinates of the topography grid (meters) that has been updated
    with the new compacted flow thickness.

- topo_gridy_initial: Vector{Float64}
    The initial y-coordinates of the topography grid (meters) prior to
    flow extrusion.

- sediment_and_flow_thickness_initial: Vector{Float64}
    The initial sediment and flow thickness (meters) prior to new flow
    extrusion.

- lava_flow_decompaction_parameters: SedimentTransportParameters
    The lava flow decompaction parameters.

- sediment_and_flow_thickness_total: Vector{Float64}
    Initialized total thickness of sediment and decompacted lava on the
    topography grid (meters).

# Updated Arrays
- topo_gridy: Vector{Float64}
    The y-coordinates of the topography grid (meters). This array is
    updated with the new de-compacted flow thickness and compacted
    pre-existing sediments and flows.

- model.markers.arrays.location.marker_y,array: Vector{Float64}
        The y-coordinates of the markers (meters). This array is updated
        using compaction field associated with flow extrusion. This is
        update occurs in the function
        apply_compaction_correction_for_topography_and_markers.

# Returns
- total_lava_thickness: Vector{Float64}
    The thickness of the new flows after decompaction (meters).

- sediment_and_flow_thickness_total: Vector{Float64}
    The total thickness of sediment and decompacted lava on the topography
    grid (meters). This includes compacted pre-existing sediment and flows
    and decompacted new flows.

- sediment_and_flow_thickness_initial_compacted: Vector{Float64}
    The initial sediment and flow thickness after compaction (meters).

"""
function apply_compaction_model(
    model::ModelData,
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_gridy_initial::Vector{Float64},
    sediment_and_flow_thickness_initial::Vector{Float64},
    lava_flow_decompaction_parameters::SedimentTransportParameters
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    (
        markers_indices_sedimentary_basin,
        markers_indices_sticky
    ) = calculate_swarm_indices_for_sediment_and_sticky(model)

    x_sorted_marker_indices_sticky = calculate_x_sorted_swarm_indices(
        model, markers_indices_sticky)

    x_sorted_marker_indices_sedimentary_basin = calculate_x_sorted_swarm_indices(
        model, markers_indices_sedimentary_basin)

    (
        total_lava_thickness,
        sediment_and_flow_thickness_total,
        sediment_and_flow_thickness_initial_compacted
    ) = apply_compaction_correction_variable_property(
        model, topo_gridx, topo_gridy, topo_gridy_initial,
        sediment_and_flow_thickness_initial,
        lava_flow_decompaction_parameters,
        markers_indices_sedimentary_basin,
        markers_indices_sticky,
        x_sorted_marker_indices_sedimentary_basin,
        x_sorted_marker_indices_sticky
    )

    return (
        total_lava_thickness,
        sediment_and_flow_thickness_total,
        sediment_and_flow_thickness_initial_compacted
    )
end

""" Decompact flow and compact pre-existing sediments and flows.

# Arguments
- model: ModelData
    The model data object.
- topo_gridx: Vector{Float64}
    The x-coordinates of the topography grid (meters).
- topo_gridy: Vector{Float64}
    The y-coordinates of the topography grid (meters) that has been updated
    with the new compacted flow thickness.
- topo_gridy_initial: Vector{Float64}
    The initial y-coordinates of the topography grid (meters) prior to
    flow extrusion.
- sediment_and_flow_thickness_initial: Vector{Float64}
    The initial sediment and flow thickness (meters) prior to new flow
    extrusion.
- sediment_transport_parameters: SedimentTransportParameters
    The sediment transport parameters.
- markers_indices_sedimentary_basin: Vector{Int64}
    The indices of the markers in the sedimentary basin.
- markers_indices_sticky: Vector{Int64}
    The indices of the sticky markers.

# Updated Arrays
- topo_gridy: Vector{Float64}
    The y-coordinates of the topography grid (meters). This array is
    updated with the new de-compacted flow thickness and compacted
    pre-existing sediments and flows.
- model.markers.arrays.location.marker_y.array: Vector{Float64}
        The y-coordinates of the markers (meters). This array is updated
        using compaction field associated with flow extrusion. This is
        update occurs in the function
        apply_compaction_correction_for_topography_and_markers.

# Returns
- new_flow_thickness_decompacted: Vector{Float64}
    The thickness of the new flow after decompaction (meters).
- sediment_and_flow_thickness_total: Vector{Float64}
    The total thickness of sediment and decompacted lava on the topography
    grid (meters). This includes compacted pre-existing sediment and flows
    and decompacted new flows.
- sediment_and_flow_thickness_initial_compacted: Vector{Float64}
    The initial sediment and flow thickness after compaction (meters).

"""
function apply_compaction_correction_variable_property(
    model::Union{ModelData, Nothing},
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_gridy_initial::Vector{Float64},
    sediment_and_flow_thickness_initial::Union{Vector{Float64}, Nothing},
    sediment_transport_parameters::SedimentTransportParameters,
    markers_indices_sedimentary_basin::Vector{Int64},
    markers_indices_sticky::Vector{Int64},
    x_sorted_marker_indices_sedimentary_basin::Vector{Int64},
    x_sorted_marker_indices_sticky::Vector{Int64}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    (
        topo_gridy_corrected,
        sediment_and_flow_thickness_total,
        new_flow_thickness_decompacted
    ) = apply_compaction_correction_for_topography_and_markers(
        model,
        sediment_transport_parameters,
        topo_gridx,
        topo_gridy_initial,
        topo_gridy,
        sediment_and_flow_thickness_initial,
        markers_indices_sedimentary_basin,
        markers_indices_sticky,
        x_sorted_marker_indices_sedimentary_basin,
        x_sorted_marker_indices_sticky
    )

    copy_array_1d!(topo_gridy_corrected, topo_gridy)

    (
        sediment_and_flow_thickness_initial_compacted
    ) = calculate_sediment_and_flow_thickness_initial_compacted(
            model, topo_gridx, markers_indices_sedimentary_basin,
            x_sorted_marker_indices_sedimentary_basin
        )

    return (
        new_flow_thickness_decompacted,
        sediment_and_flow_thickness_total,
        sediment_and_flow_thickness_initial_compacted
    )
end

function calculate_sediment_and_flow_thickness_initial_compacted(
    model::ModelData,
    topo_gridx::Vector{Float64},
    markers_indices_sedimentary_basin::Vector{Int64},
    x_sorted_marker_indices_sedimentary_basin::Vector{Int64}
)::Vector{Float64}
    search_radius = topo_gridx[2] - topo_gridx[1]
    markers_x = model.markers.arrays.location.marker_x.array
    markers_y = model.markers.arrays.location.marker_y.array

    nsmooth_top_bottom = model.topography.parameters.topo_grid.nsmooth_top_bottom.value

    (
        top_marker_sediment, 
        bottom_marker_sediment
    ) = calculate_top_and_bottom_of_swarm_opt(
        x_sorted_marker_indices_sedimentary_basin,
        markers_x, markers_y, topo_gridx,
        search_radius; use_smoothing=true, nsmooth=nsmooth_top_bottom
    )

    (
        sediment_and_flow_thickness_initial_compacted
    ) = bottom_marker_sediment - top_marker_sediment
    return sediment_and_flow_thickness_initial_compacted
end

end # module