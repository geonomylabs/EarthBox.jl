"""
Module for marker compaction scratch buffers.

Holds pre-allocated marker-sized buffers used by sediment-transport
compaction routines (`MarkerCompaction.compact_sediment_and_advect_markers!`
and `CompactionCorrection.apply_compaction_correction_for_topography_and_markers`).
Each buffer is sized to `marknum`; values are written/overwritten per call.
Consumers must bound their reads by the marker-index arrays passed alongside
(e.g., `markers_indices_sedimentary_basin`), not by `length(buffer)`.
"""
module CompactionGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "compaction"

const ADATA = get_eb_arrays()

"""
    Compaction <: AbstractArrayGroup

Array group holding scratch buffers reused across calls to sediment-
transport compaction routines. Replaces per-call `Vector{Int64}(undef, marknum)`
and `zeros(Float64, nmarkers)` allocations.

# Fields
- `markers_topo_xindex::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.markers_topo_xindex.description)
- `markers_compaction_yindex::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.markers_compaction_yindex.description)
- `markers_unit_distance_from_cell_top::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.markers_unit_distance_from_cell_top.description)
- `total_marker_compaction_displacement::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.total_marker_compaction_displacement.description)

# Nested Dot Access
- `model.markers.arrays.compaction.markers_topo_xindex.array`
- `model.markers.arrays.compaction.markers_compaction_yindex.array`
- `model.markers.arrays.compaction.markers_unit_distance_from_cell_top.array`
- `model.markers.arrays.compaction.total_marker_compaction_displacement.array`

# Constructor
    Compaction(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (all buffers sized to this length).
"""
mutable struct Compaction <: AbstractArrayGroup
    markers_topo_xindex::MarkerArrayInt1DState{Int64}
    markers_compaction_yindex::MarkerArrayInt1DState{Int64}
    markers_unit_distance_from_cell_top::MarkerArrayFloat1DState{Float64}
    total_marker_compaction_displacement::MarkerArrayFloat1DState{Float64}
    sticky_displacement_factors::MarkerArrayFloat1DState{Float64}
    sticky_marker_displacement::MarkerArrayFloat1DState{Float64}
end

function Compaction(marknum::Int)::Compaction
    return Compaction(
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.markers_topo_xindex.name,
            ADATA.markers_topo_xindex.units,
            ADATA.markers_topo_xindex.description
        ),
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.markers_compaction_yindex.name,
            ADATA.markers_compaction_yindex.units,
            ADATA.markers_compaction_yindex.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.markers_unit_distance_from_cell_top.name,
            ADATA.markers_unit_distance_from_cell_top.units,
            ADATA.markers_unit_distance_from_cell_top.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.total_marker_compaction_displacement.name,
            ADATA.total_marker_compaction_displacement.units,
            ADATA.total_marker_compaction_displacement.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.sticky_displacement_factors.name,
            ADATA.sticky_displacement_factors.units,
            ADATA.sticky_displacement_factors.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.sticky_marker_displacement.name,
            ADATA.sticky_marker_displacement.units,
            ADATA.sticky_marker_displacement.description
        )
    )
end

end # module
