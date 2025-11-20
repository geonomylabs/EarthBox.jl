module TopBoundary

import EarthBox.ModelDataContainer: ModelData
import ..InflowOutflowBC
import ..FreeSlipBC
import ..NoSlipBC
import ..TemperatureBC
import ..ConstantGradientBC
import ..SetBCTools

function free_slip!(model::ModelData)::Nothing
    FreeSlipBC.set_free_slip_bc_at_top_boundary!(model)
    return nothing
end

function no_slip!(model::ModelData)::Nothing
    NoSlipBC.set_no_slip_bc_at_top_boundary!(model)
    return nothing
end

function inflow_outflow!(model::ModelData, vy::Float64)::Nothing
    InflowOutflowBC.set_inflow_outflow_bc_at_top_boundary!(model, vy)
    return nothing
end

function temperature!(
    model::ModelData, 
    t_top::Union{Float64, Nothing}=nothing
)::Nothing
    if t_top === nothing
        t_top = SetBCTools.get_top_temperature(model)
    end
    TemperatureBC.set_temperature_at_top_boundary!(t_top, model)
    return nothing
end

""" 
    constant_vertical_temperature_gradient!(
        model::ModelData; 
        ystp::Union{Float64, Nothing}=nothing, 
        dtdy::Union{Float64, Nothing}=nothing
    )::Nothing

Set constant vertical temperature gradient at top boundary.

Temperature is assumed to increase with depth starting at basic nodes
located below the top boundary. Temperature increases linearly with
depth according to the defined gradient dtdy. The vertical spacing ystp
is the thickness of the uniform basic cell and is used to calculate the
temperature change across the cell via delta_temperature = dtdy*ystp.

"""
function constant_vertical_temperature_gradient!(
    model::ModelData;
    ystp::Union{Float64, Nothing}=nothing,
    dtdy::Union{Float64, Nothing}=nothing
)::Nothing
    ystp = SetBCTools.check_temperature_gradient_and_update_spacing(dtdy, ystp)
    if ystp === nothing || dtdy === nothing
        ystp, dtdy = SetBCTools.calculate_vertical_gradient(model)
    end

    ConstantGradientBC.set_constant_vertical_gradient_at_top_boundary!(
        ystp, dtdy, model)
    return nothing
end

end # module