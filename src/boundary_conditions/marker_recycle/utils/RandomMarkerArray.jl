module RandomMarkerArray

import Random: rand!
import EarthBox.ModelDataContainer: ModelData

""" Refill and return the persistent marknum-sized random-number buffer.

The buffer (`model.markers.arrays.solidification.marker_random_buffer.array`)
is shared across pre-solver `Solidification.solidify!` and post-solver
`MarkerRecycle.ResetMarkers.reset_subsurface_markers!`. Each consumer is
responsible for calling this function (or directly `rand!` on the buffer)
to refresh the values before reading. The two consumers run sequentially
within a single time step, so values do not leak between them.

The returned reference points to the persistent buffer; callers must treat
it as read-only and not store the reference past their own scope.
"""
function get_random_marker_array(model::ModelData)::Vector{Float64}
    buf = model.markers.arrays.solidification.marker_random_buffer.array
    rand!(buf)
    return buf
end

end # module