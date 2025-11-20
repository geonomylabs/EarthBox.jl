module LeftBoundary

import EarthBox.ModelDataContainer: ModelData
import ..InflowOutflowBC
import ..FreeSlipBC
import ..NoSlipBC
import ..TemperatureBC
import ..ZeroHeatFluxBC
import ..SetBCTools

function free_slip!(model::ModelData)::Nothing
    FreeSlipBC.set_free_slip_bc_at_left_boundary!(model)
    return nothing
end

function no_slip!(model::ModelData)::Nothing
    NoSlipBC.set_no_slip_bc_at_left_boundary!(model)
    return nothing
end

function inflow_outflow!(model::ModelData, vx::Float64)::Nothing
    InflowOutflowBC.set_inflow_outflow_bc_at_left_boundary!(model, vx)
    return nothing
end

function inflow_outflow_depth_dependent!(
    model::ModelData, 
    y_linear::Float64, 
    vx::Float64
)::Nothing
    InflowOutflowBC.set_inflow_outflow_bc_at_left_boundary_depth_dependent!(
        model, y_linear, vx)
    return nothing
end

function inflow_and_outflow_along_side!(
    model::ModelData,
    plate_velocity::Float64,
    inflow_velocity::Float64
)::Nothing
    InflowOutflowBC.set_inflow_and_outflow_bc_at_left_boundary!(
        model, plate_velocity, inflow_velocity)
    return nothing
end

function temperature!(
    model::ModelData;
    t_left::Union{Float64, Nothing}=nothing
)::Nothing
    if t_left === nothing
        t_left = SetBCTools.get_left_temperature(model)
    end
    TemperatureBC.set_temperature_at_left_boundary!(t_left, model)
    return nothing
end

""" 
    variable_temperature!(
        model::ModelData; 
        t_top::Union{Float64, Nothing}=nothing, 
        ystp::Union{Float64, Nothing}=nothing, 
        dtdy::Union{Float64, Nothing}=nothing
    )::Nothing

Temperature increases linearly with depth starting at the top of the
model where temperature is t_top and changing according to the defined
gradient dtdy. The vertical spacing ystp is used to calculate depth at
each vertical grid point via float(i)*ystp where i is the index of the
grid node assuming uniform grid spacing.

"""
function variable_temperature!(
    model::ModelData;
    t_top::Union{Float64, Nothing}=nothing,
    ystp::Union{Float64, Nothing}=nothing,
    dtdy::Union{Float64, Nothing}=nothing
)::Nothing
    if ystp === nothing || dtdy === nothing
        ystp, dtdy = SetBCTools.calculate_vertical_gradient(model)
    end
    if t_top === nothing
        t_top = SetBCTools.get_top_temperature(model)
    end
    SetBCTools.check_for_uniform_grid_spacing(
        model, "set_variable_temperature_at_left_boundary")
    TemperatureBC.set_variable_temperature_at_left_boundary!(
        t_top, ystp, dtdy, model)
    return nothing
end

function zero_heat_flux!(model::ModelData)::Nothing
    ZeroHeatFluxBC.set_zero_heat_flux_at_left_boundary!(model)
    return nothing
end

end # module LeftBoundary 