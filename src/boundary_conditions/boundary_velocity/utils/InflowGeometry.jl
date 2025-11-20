module InflowGeometry

import EarthBox.ModelDataContainer: ModelData
import ..OutflowGeometry: get_outflow_geometry
import ..StickyThickness: get_sticky_thickness

"""
    get_y_height_of_inflow(model::ModelData, y_initial::Float64)::Float64

Get height of inflow zone along side boundaries.

# Arguments
- `model::ModelData`: Model data container
- `y_initial::Float64`: Initial y-coordinate

# Returns
- Height of inflow zone
"""
function get_y_height_of_inflow(
    model::ModelData,
    y_initial::Float64
)::Float64
    ysize = model.grids.parameters.geometry.ysize.value
    return ysize - y_initial
end

"""
    get_starting_inflow_coordinates_left(model::ModelData)::Tuple{Float64, Float64}

Get starting coordinates for recycling markers on left side.

# Arguments
- `model::ModelData`: Model data container

# Returns
- Tuple containing:
  - `x_left_initial::Float64`: Initial x-coordinate for left side inflow boundary
  - `y_left_initial::Float64`: Initial y-coordinate for top of left side inflow boundary
"""
function get_starting_inflow_coordinates_left(
    model::ModelData
)::Tuple{Float64, Float64}
    sticky_thickness_left, _ = get_sticky_thickness(model)
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    x_left_initial = 0.0
    y_left_initial = calculate_top_of_side_inflow(
        sticky_thickness_left, plate_thickness, smoothing_thickness)
    return (x_left_initial, y_left_initial)
end

"""
    get_starting_inflow_coordinates_right(model::ModelData)::Tuple{Float64, Float64}

Get starting coordinates for recycling markers on right side.

# Arguments
- `model::ModelData`: Model data container

# Returns
- Tuple containing:
  - `x_right_initial::Float64`: Initial x-coordinate for right side inflow boundary
  - `y_right_initial::Float64`: Initial y-coordinate for top right side inflow boundary
"""
function get_starting_inflow_coordinates_right(
    model::ModelData
)::Tuple{Float64, Float64}
    _, sticky_thickness_right = get_sticky_thickness(model)
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    x_right_initial = model.grids.parameters.geometry.xsize.value
    y_right_initial = calculate_top_of_side_inflow(
        sticky_thickness_right, plate_thickness, smoothing_thickness)
    return (x_right_initial, y_right_initial)
end

"""
    calculate_top_of_side_inflow(sticky_thickness::Float64, plate_thickness::Float64, smoothing_thickness::Float64)::Float64

Calculate y-depth of the top of the side inflow zone.

# Arguments
- `sticky_thickness::Float64`: Thickness of sticky layer
- `plate_thickness::Float64`: Thickness of plate
- `smoothing_thickness::Float64`: Thickness of smoothing layer

# Returns
- `::Float64`: Depth of the top of the side inflow zone
"""
function calculate_top_of_side_inflow(
    sticky_thickness::Float64,
    plate_thickness::Float64,
    smoothing_thickness::Float64
)::Float64
    return sticky_thickness + plate_thickness + smoothing_thickness
end

"""
    get_inflow_geometry(model::ModelData)::Tuple{Float64, Float64}

Get inflow geometry.

# Arguments
- `model::ModelData`: Model data container

# Returns
- Tuple containing:
  - `ysize_inflow_left::Float64`: Height (m) of inflow zone on left side
  - `ysize_inflow_right::Float64`: Height (m) of inflow zone on right side
"""
function get_inflow_geometry(
    model::ModelData
)::Tuple{Float64, Float64}
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    sticky_thickness_left, sticky_thickness_right = get_sticky_thickness(model)
    ysize = model.grids.parameters.geometry.ysize.value
    ysize_inflow_left = (
        ysize - sticky_thickness_left - plate_thickness - smoothing_thickness
    )
    ysize_inflow_right = (
        ysize - sticky_thickness_right - plate_thickness - smoothing_thickness
    )
    return (ysize_inflow_left, ysize_inflow_right)
end

end # module InflowGeometry 