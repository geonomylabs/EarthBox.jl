module InflowOutflowBC

using EarthBox.ModelDataContainer: ModelData
using ..StokesGradientBC
using ..VelocityBC
using ..VelocityInflowOutflow

function set_inflow_outflow_bc_at_top_boundary!(
    model::ModelData,
    vy::Float64
)::Nothing
    StokesGradientBC.set_zero_vertical_vx_gradient_at_top_boundary!(model)
    VelocityBC.set_velocity_y_at_top_boundary_ghost_nodes!(vy, model)
    return nothing
end

function set_inflow_outflow_bc_at_bottom_boundary!(
    model::ModelData,
    vy::Float64
)::Nothing
    StokesGradientBC.set_zero_vertical_vx_gradient_at_bottom_boundary!(model)
    VelocityBC.set_velocity_y_at_bottom_boundary_ghost_nodes!(vy, model)
    return nothing
end

function set_inflow_outflow_bc_at_left_boundary!(
    model::ModelData,
    vx::Float64
)::Nothing
    StokesGradientBC.set_zero_lateral_vy_gradient_at_left_boundary!(model)
    VelocityBC.set_velocity_x_at_left_boundary_ghost_nodes!(vx, model)
    return nothing
end

function set_inflow_outflow_bc_at_right_boundary!(
    model::ModelData,
    vx::Float64
)::Nothing
    StokesGradientBC.set_zero_lateral_vy_gradient_at_right_boundary!(model)
    VelocityBC.set_velocity_x_at_right_boundary_ghost_nodes!(vx, model)
    return nothing
end

function set_inflow_outflow_bc_at_left_boundary_depth_dependent!(
    model::ModelData,
    y_linear::Float64,
    vx::Float64
)::Nothing
    bc_gridy = get_boundary_condition_vertical_grid_nodes(model)
    StokesGradientBC.set_zero_lateral_vy_gradient_at_left_boundary!(model)
    VelocityBC.set_linear_velocity_x_at_left_boundary_ghost_nodes!(
        y_linear, vx, bc_gridy, model)
    return nothing
end

function set_inflow_outflow_bc_at_right_boundary_depth_dependent!(
    model::ModelData,
    y_linear::Float64,
    vx::Float64
)::Nothing
    bc_gridy = get_boundary_condition_vertical_grid_nodes(model)
    StokesGradientBC.set_zero_lateral_vy_gradient_at_right_boundary!(model)
    VelocityBC.set_linear_velocity_x_at_right_boundary_ghost_nodes!(
        y_linear, vx, bc_gridy, model)
    return nothing
end

function set_inflow_and_outflow_bc_at_left_boundary!(
    model::ModelData,
    plate_velocity::Float64,
    inflow_velocity::Float64
)::Nothing
    bc_gridy = get_boundary_condition_vertical_grid_nodes(model)
    StokesGradientBC.set_zero_lateral_vy_gradient_at_left_boundary!(model)
    VelocityInflowOutflow.set_inflow_and_outflow_velocity_x_at_left_boundary_ghost_nodes!(
        plate_velocity, inflow_velocity, bc_gridy, model)
    return nothing
end

function set_inflow_and_outflow_bc_at_right_boundary!(
    model::ModelData,
    plate_velocity::Float64,
    inflow_velocity::Float64
)::Nothing
    bc_gridy = get_boundary_condition_vertical_grid_nodes(model)
    StokesGradientBC.set_zero_lateral_vy_gradient_at_right_boundary!(model)
    VelocityInflowOutflow.set_inflow_and_outflow_velocity_x_at_right_boundary_ghost_nodes!(
        plate_velocity, inflow_velocity, bc_gridy, model)
    return nothing
end

function get_boundary_condition_vertical_grid_nodes(
    model::ModelData
)::Vector{Float64}
    ynum = model.grids.parameters.geometry.ynum.value
    gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array
    bc_gridy = zeros(Float64, ynum+1)
    for i in 1:ynum+1
        bc_gridy[i] = gridy_vx[i]
    end
    return bc_gridy
end

end # module