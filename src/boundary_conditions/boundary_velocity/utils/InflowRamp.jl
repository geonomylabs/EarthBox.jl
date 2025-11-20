module InflowRamp

import EarthBox.ModelDataContainer: ModelData
import ..StickyThickness: get_sticky_thickness
import ..OutflowGeometry: get_outflow_geometry

"""
    get_top_and_bottom_of_left_or_right_ramp(
        model::ModelData; left_side::Bool=false, right_side::Bool=false
    )::Tuple{Float64, Float64}

Get top and bottom of left or right inflow ramp.

# Arguments
- `model::ModelData`: Model data container
- `left_side::Bool=false`: If true, calculate for left side
- `right_side::Bool=false`: If true, calculate for right side

# Returns
- Tuple containing:
  - y_smooth_top: Top of inflow ramp
  - y_smooth_bottom: Bottom of inflow ramp
"""
function get_top_and_bottom_of_left_or_right_ramp(
    model::ModelData;
    left_side::Bool=false,
    right_side::Bool=false
)::Tuple{Float64, Float64}
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    sticky_thickness_left, sticky_thickness_right = get_sticky_thickness(model)
    if left_side
        sticky_thickness = sticky_thickness_left
    elseif right_side
        sticky_thickness = sticky_thickness_right
    else
        sticky_thickness = sticky_thickness_left
        println(
            "!!! WARNING !!!: No side specified. Defaulting to left side " *
            "inflow ramp."
        )
    end
    y_smooth_top, y_smooth_bottom = get_top_and_bottom_of_inflow_ramp(
        sticky_thickness, plate_thickness, smoothing_thickness)
    return (y_smooth_top, y_smooth_bottom)
end

"""
    calculate_ramp_factor(
        y_marker::Float64, y_top_ramp::Float64, y_bottom_ramp::Float64)::Float64

Calculate ramp factor used to smooth velocity transition.

The transition from outflow to inflow is smoothed out using a ramp.

# Arguments
- `y_marker::Float64`: Y-coordinate of marker
- `y_top_ramp::Float64`: Top of ramp
- `y_bottom_ramp::Float64`: Bottom of ramp

# Returns
- Linear factor for ramp based on depth:
  - 0.0 = top of ramp
  - 0.5 = middle of ramp
  - 1.0 = base of ramp
  - If marker is located outside of ramp, ramp_factor = 1.0
"""
function calculate_ramp_factor(
    y_marker::Float64,
    y_top_ramp::Float64,
    y_bottom_ramp::Float64
)::Float64
    if y_top_ramp <= y_marker < y_bottom_ramp
        ramp_factor = (y_marker - y_top_ramp)/(y_bottom_ramp - y_top_ramp)
    else
        ramp_factor = 1.0
    end
    return ramp_factor
end

"""
    get_top_and_bottom_of_inflow_ramp(
        sticky_thickness::Float64, 
        plate_thickness::Float64, 
        smoothing_thickness::Float64
    )::Tuple{Float64, Float64}

Get top and bottom of inflow ramp.

# Arguments
- `sticky_thickness::Float64`: Thickness of sticky layer
- `plate_thickness::Float64`: Thickness of plate
- `smoothing_thickness::Float64`: Thickness of smoothing layer

# Returns
- Tuple containing top and bottom of inflow ramp
"""
function get_top_and_bottom_of_inflow_ramp(
    sticky_thickness::Float64,
    plate_thickness::Float64,
    smoothing_thickness::Float64
)::Tuple{Float64, Float64}
    y_smooth_top = sticky_thickness + plate_thickness + smoothing_thickness
    y_smooth_bottom = y_smooth_top + smoothing_thickness
    return (y_smooth_top, y_smooth_bottom)
end

end # module