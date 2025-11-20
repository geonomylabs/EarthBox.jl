module OutflowGeometry

import EarthBox.ModelDataContainer: ModelData
import ..StickyThickness: get_sticky_thickness

function get_outflow_geometry(
    model::ModelData
)::Tuple{Float64, Float64}
    plate_thickness = model.bcs.parameters.velocity.plate_thickness.value
    smoothing_thickness = model.bcs.parameters.velocity.smoothing_thickness.value
    return (plate_thickness, smoothing_thickness)
end

function get_outflow_size(
    model::ModelData
)::Tuple{Float64, Float64}
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    sticky_thickness_left, sticky_thickness_right = get_sticky_thickness(model)

    outflow_size_right = (
        sticky_thickness_right + plate_thickness + smoothing_thickness)

    outflow_size_left = (
        sticky_thickness_left + plate_thickness + smoothing_thickness)

    return (outflow_size_left, outflow_size_right)
end

end # module OutflowGeometry 