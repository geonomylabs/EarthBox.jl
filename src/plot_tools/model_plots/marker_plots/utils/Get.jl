module Get

import ...PlotParametersManager: PlotParameters
import ..PlotMarkerArraysManager: PlotMarkerArrays

function get_color_bar_axis_fractions(
    parameters::PlotParameters
)::Tuple{Float64, Float64}
    axis_fraction_for_color_bar = parameters.color_bar.axis_fraction_for_color_bar
    axis_fraction_for_color_bar_gap = parameters.color_bar.axis_fraction_for_color_bar_gap
    return axis_fraction_for_color_bar, axis_fraction_for_color_bar_gap
end

function get_marker_coordinates(
    marker_arrays::PlotMarkerArrays
)::Tuple{Vector{Float64}, Vector{Float64}}
    x_array_km = copy(marker_arrays.marker_x_km)
    y_array_km = copy(marker_arrays.marker_y_km)
    return x_array_km, y_array_km
end

end # module
