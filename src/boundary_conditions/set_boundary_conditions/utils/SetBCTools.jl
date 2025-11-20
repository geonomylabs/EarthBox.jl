module SetBCTools

import EarthBox.ModelDataContainer: ModelData

function calculate_vertical_gradient(model::ModelData)::Tuple{Float64, Float64}
    check_for_uniform_grid_spacing(model, "calculate_vertical_gradient")
    ynum = model.grids.parameters.geometry.ynum.value
    ysize = model.grids.parameters.geometry.ysize.value

    temperature_parameters = model.bcs.parameters.temperature
    t_top = temperature_parameters.temperature_top.value
    t_bottom = temperature_parameters.temperature_bottom.value

    ystp = ysize/float(ynum-1)
    dtdy = (t_bottom - t_top)/ysize

    return ystp, dtdy
end

function get_top_temperature(model::ModelData)::Float64
    temperature_parameters = model.bcs.parameters.temperature
    return temperature_parameters.temperature_top.value
end

function get_bottom_temperature(model::ModelData)::Float64
    temperature_parameters = model.bcs.parameters.temperature
    return temperature_parameters.temperature_bottom.value
end

function get_right_temperature(model::ModelData)::Float64
    temperature_parameters = model.bcs.parameters.temperature
    return temperature_parameters.temperature_right.value
end

function get_left_temperature(model::ModelData)::Float64
    temperature_parameters = model.bcs.parameters.temperature
    return temperature_parameters.temperature_left.value
end

function check_temperature_gradient_and_update_spacing(
    dtdy::Union{Float64, Nothing}, 
    ystp::Union{Float64, Nothing}
)::Union{Float64, Nothing}
    if dtdy == 0.0
        ystp = 1.0
    end
    return ystp
end

function check_for_uniform_grid_spacing(
    model::ModelData, 
    description::String
)::Nothing
    itype_grid = model.grids.parameters.grid_options.itype_grid.value
    if itype_grid != 0
        stype_grid = model.grids.parameters.grid_options.stype_grid.value
        error(
            "$description is only implemented for uniform grid spacing (itype_grid = 0). "
            *"The current grid type is itype_grid = $itype_grid ($stype_grid)."
            )
    end
    return nothing
end

end # module SetBCTools 