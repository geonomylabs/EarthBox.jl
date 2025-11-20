module PartialMeltIndices

import EarthBox.ModelDataContainer: ModelData

""" Get indices of partially molten markers.

Partial melting is defined as markers with melt fraction greater than 0.0
and extracted melt fraction greater than or equal to 0.0. This condition
is assumed to indicate that the marker is partially molten and channelized
porosity is available for melt migration.

# Returns
- partial_melt_marker_indices::Vector{Int64}
    - Contains indices for all partially molten particles.
- nmarkers_partial_melt::Int
    - Number of partially molten markers.
"""
function get_partial_melt_marker_indices(
    model::ModelData;
    xstart::Float64=-1e39,
    xend::Float64=1e39
)::Tuple{Vector{Int64}, Int}
    marknum = model.markers.parameters.distribution.marknum.value
    
    partial_melt_marker_indices = Vector{Int64}(undef, marknum) #zeros(Int64, marknum)
    Threads.@threads for imarker in 1:marknum
        partial_melt_marker_indices[imarker] = 0
    end

    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
    marker_extracted_meltfrac = model.markers.arrays.melt.marker_extracted_meltfrac.array
    marker_x = model.markers.arrays.location.marker_x.array

    icount = 0
    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        if xstart <= x_marker <= xend
            meltfrac = marker_meltfrac[imarker]
            extracted_meltfrac = marker_extracted_meltfrac[imarker]
            if meltfrac > 0.0 && extracted_meltfrac >= 0.0
                partial_melt_marker_indices[icount + 1] = imarker
                icount += 1
            end
        end
    end
    nmarkers_partial_melt = icount
    return partial_melt_marker_indices, nmarkers_partial_melt
end

end # module 