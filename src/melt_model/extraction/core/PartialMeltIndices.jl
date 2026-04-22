module PartialMeltIndices

import EarthBox.ModelDataContainer: ModelData

""" Update the packed list of partially molten marker indices.

Partial melting is defined as markers with melt fraction greater than 0.0
and extracted melt fraction greater than or equal to 0.0. This condition
is assumed to indicate that the marker is partially molten and channelized
porosity is available for melt migration.

Writes the packed prefix `[1:nmarkers_partial_melt]` into the pre-allocated
scratch buffer at `model.melting.arrays.buffers.partial_melt_marker_indices.array`.
Positions past `nmarkers_partial_melt` carry stale values from prior calls and
must not be read.

# Returns
- nmarkers_partial_melt::Int
    - Number of partially molten markers.
"""
function update_partial_melt_marker_indices!(
    model::ModelData;
    xstart::Float64=-1e39,
    xend::Float64=1e39
)::Int
    marknum = model.markers.parameters.distribution.marknum.value
    partial_melt_marker_indices =
        model.melting.arrays.buffers.partial_melt_marker_indices.array

    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
    marker_extracted_meltfrac = model.markers.arrays.melt.marker_extracted_meltfrac.array
    marker_x = model.markers.arrays.location.marker_x.array

    icount = 0
    @inbounds for imarker in 1:marknum
        x_marker = marker_x[imarker]
        if xstart <= x_marker <= xend
            meltfrac = marker_meltfrac[imarker]
            extracted_meltfrac = marker_extracted_meltfrac[imarker]
            if meltfrac > 0.0 && extracted_meltfrac >= 0.0
                icount += 1
                partial_melt_marker_indices[icount] = imarker
            end
        end
    end
    return icount
end

end # module
