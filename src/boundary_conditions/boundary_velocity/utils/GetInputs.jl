module GetInputs

import EarthBox.ModelDataContainer: ModelData
import ..StickyThickness: get_sticky_thickness

struct ExtensionBCInputs
    velocity_left::Float64
    velocity_right::Float64
    sticky_thickness_left::Float64
    sticky_thickness_right::Float64
    xsize::Float64
    ysize::Float64
end

function get_extension_inputs(
    model::ModelData;
    use_asymmetric_inflow_outflow::Bool=false
)::ExtensionBCInputs
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value

    sticky_thickness_left, sticky_thickness_right = get_sticky_thickness(model)

    velocity_left, velocity_right = get_left_and_right_extension_velocities(
        model, use_asymmetric_inflow_outflow=use_asymmetric_inflow_outflow)

    return ExtensionBCInputs(
        velocity_left, velocity_right,
        sticky_thickness_left, sticky_thickness_right,
        xsize, ysize
    )
end

function get_left_and_right_extension_velocities(
    model::ModelData;
    use_asymmetric_inflow_outflow::Bool=false
)::Tuple{Float64, Float64}
    full_velocity_extension = abs(model.bcs.parameters.velocity.full_velocity_extension.value)
    if !use_asymmetric_inflow_outflow
        velocity_left = -full_velocity_extension/2.0
        velocity_right = full_velocity_extension/2.0
    else
        velocity_left = 0.0
        velocity_right = full_velocity_extension
    end
    return (velocity_left, velocity_right)
end

end # module