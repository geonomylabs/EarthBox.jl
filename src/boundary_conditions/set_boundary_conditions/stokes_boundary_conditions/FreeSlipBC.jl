module FreeSlipBC

using EarthBox.ModelDataContainer: ModelData
using ..StokesGradientBC
using ..VelocityBC

function set_free_slip_bc_at_top_boundary!(model::ModelData)::Nothing
    StokesGradientBC.set_zero_vertical_vx_gradient_at_top_boundary!(model)
    VelocityBC.set_velocity_y_at_top_boundary_ghost_nodes!(0.0, model)
    return nothing
end

function set_free_slip_bc_at_bottom_boundary!(model::ModelData)::Nothing
    StokesGradientBC.set_zero_vertical_vx_gradient_at_bottom_boundary!(model)
    VelocityBC.set_velocity_y_at_bottom_boundary_ghost_nodes!(0.0, model)
    return nothing
end

function set_free_slip_bc_at_left_boundary!(model::ModelData)::Nothing
    StokesGradientBC.set_zero_lateral_vy_gradient_at_left_boundary!(model)
    VelocityBC.set_velocity_x_at_left_boundary_ghost_nodes!(0.0, model)
    return nothing
end

function set_free_slip_bc_at_right_boundary!(model::ModelData)::Nothing
    StokesGradientBC.set_zero_lateral_vy_gradient_at_right_boundary!(model)
    VelocityBC.set_velocity_x_at_right_boundary_ghost_nodes!(0.0, model)
    return nothing
end

end # module FreeSlipBC 