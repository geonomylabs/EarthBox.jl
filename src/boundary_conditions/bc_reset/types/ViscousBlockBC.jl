module ViscousBlockBC

import EarthBox.ModelDataContainer: ModelData
import ...SetBoundaryConditions: StokesGradientBC
import ...SetBoundaryConditions: VelocityBC
import ...SetBoundaryConditions: PressureBC
import ...SetBoundaryConditions: TemperatureBC

function set_bcs!(model::ModelData)
    # Set pressure boundary condition to upper right corner mode
    PressureBC.set_pressure_bc_mode!(model, 0)
    # Velocity boundary conditions along top boundary
    StokesGradientBC.set_zero_vertical_vx_gradient_at_top_boundary!(model)
    VelocityBC.set_velocity_y_at_top_boundary_ghost_nodes!(0.0, model)
    # Velocity boundary condition along bottom boundary
    StokesGradientBC.set_zero_vertical_vx_gradient_at_bottom_boundary!(model)
    VelocityBC.set_velocity_y_at_bottom_boundary_ghost_nodes!(0.0, model)
    # Velocity boundary conditions along left boundary
    VelocityBC.set_velocity_x_at_left_boundary_ghost_nodes!(0.0, model)
    StokesGradientBC.set_zero_lateral_vy_gradient_at_left_boundary!(model)
    # Velocity boundary conditions along right boundary
    VelocityBC.set_velocity_x_at_right_boundary_ghost_nodes!(0.0, model)
    StokesGradientBC.set_zero_lateral_vy_gradient_at_right_boundary!(model)
    # Temperature boundary condition along top boundary
    TemperatureBC.set_temperature_at_top_boundary!(0.0, model)
    # Temperature boundary condition along bottom boundary
    TemperatureBC.set_temperature_at_bottom_boundary!(0.0, model)
    # Temperature condition along left boundary
    TemperatureBC.set_temperature_at_left_boundary!(0.0, model)
    # Temperature boundary condition along right boundary
    TemperatureBC.set_temperature_at_right_boundary!(0.0, model)
end

end # module ViscousBlockBC 