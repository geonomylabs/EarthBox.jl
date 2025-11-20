module TimeStepUtils

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: years_to_seconds, seconds_to_years
import EarthBox.StaggeredGrid: check_for_ttype_refinement

function calculate_time_step_parameters(
    characteristic_max_velocity_m_s::Float64,
    minimum_cell_spacing_m::Float64,
    marker_cell_displ_max::Float64,
    model_duration_myr::Float64,
    buffer_time_myr::Float64,
    timestep_viscoelastic_limit_yrs::Float64
)::Tuple{Int, Float64}
    timestep_viscoelastic_yrs = calculate_viscoelastic_time_step(
        characteristic_max_velocity_m_s, minimum_cell_spacing_m, 
        marker_cell_displ_max, timestep_viscoelastic_limit_yrs
        )
    ntimestep_max = calculate_number_of_timesteps(
        timestep_viscoelastic_yrs, model_duration_myr, buffer_time_myr)
    timestep_viscoelastic_s = years_to_seconds(timestep_viscoelastic_yrs)
    return ntimestep_max, timestep_viscoelastic_s
end

function calculate_viscoelastic_time_step(
    characteristic_max_velocity_m_s::Float64,
    minimum_cell_spacing_m::Float64,
    marker_cell_displ_max::Float64,
    timestep_viscoelastic_limit_yrs::Float64
)::Float64
    marker_displacement_limit_m = minimum_cell_spacing_m * marker_cell_displ_max
    if isnan(characteristic_max_velocity_m_s)
        return NaN
    end
    (
        timestep_viscoelastic_s
    ) = marker_displacement_limit_m / characteristic_max_velocity_m_s
    timestep_viscoelastic_yrs = seconds_to_years(timestep_viscoelastic_s)
    print_info("Calculated timestep_viscoelastic_yrs: $(timestep_viscoelastic_yrs) (yrs)", level=1)
    timestep_viscoelastic_yrs = min(timestep_viscoelastic_yrs, timestep_viscoelastic_limit_yrs)
    return timestep_viscoelastic_yrs
end

function calculate_number_of_timesteps(
    timestep_viscoelastic_yrs::Float64,
    model_duration_myr::Float64,
    buffer_time_myr::Float64
)::Int
    ntimestep_max = floor(Int, model_duration_myr * 1e6 / timestep_viscoelastic_yrs)
    # Add a small time buffer to ensure sufficient time steps given variability from 
    # adaptive time stepping
    if buffer_time_myr > 0.0
        ntimestep_max += floor(Int, buffer_time_myr * 1e6 / timestep_viscoelastic_yrs)
    end
    return ntimestep_max
end

""" Calculate number of output time steps.

This applies to cases where there is no step in velocity boundary conditions
and visco-elastic time step.

If iend is defined then the number of time steps is calculated as the
difference between iend and istart plus one taking into account a one 
starting point.

# Arguments
- `timestep_viscoelastic_yrs::Float64`: Viscoelastic time step in years
- `timestep_out_yrs::Float64`: Output time step in years
- `model_duration_myr::Float64`: Model duration in Myr
- `istart::Int`: Index of the first time step (default: 1)
- `iend::Union{Int, Nothing}`: Index of the last time step (default: nothing)

# Returns
- `number_of_output_time_steps::Float64`: Number of output time steps
"""
function calculate_number_of_output_time_steps(
    timestep_viscoelastic_yrs::Float64,
    timestep_out_yrs::Float64,
    model_duration_myr::Float64,
    istart::Int=1,
    iend::Union{Int, Nothing}=nothing
)::Float64
    if iend === nothing
        model_duration_yrs = model_duration_myr * 1e6
        ntimestep_max = floor(Int, 
            model_duration_yrs / timestep_viscoelastic_yrs /
            floor(Int, timestep_out_yrs / timestep_viscoelastic_yrs))
    else
        ntimestep_max = (iend - istart) + 1
    end
    return ntimestep_max
end

function calculate_minimum_cell_spacing(model::ModelData)::Float64
    if check_for_ttype_refinement(model)
        minimum_cell_spacing_m = calculate_minimum_cell_spacing_for_ttype_refinement(
            model.grids.parameters.refinement.dx_highres.value,
            model.grids.parameters.refinement.dy_highres.value
        )
    else
        minimum_cell_spacing_m = calculate_minimum_cell_spacing_for_general_grid(
            model.grids.parameters.geometry.xsize.value,
            model.grids.parameters.geometry.ysize.value,
            model.grids.parameters.geometry.xnum.value,
            model.grids.parameters.geometry.ynum.value
        )
    end
end

function calculate_minimum_cell_spacing_for_ttype_refinement(
    dx_highres_m::Float64,
    dy_highres_m::Float64,
)::Float64
    minimum_cell_spacing_m = min(dx_highres_m, dy_highres_m)
    return minimum_cell_spacing_m
end

function calculate_minimum_cell_spacing_for_general_grid(
    xsize::Float64,
    ysize::Float64,
    xnum::Int64,
    ynum::Int64
)::Float64
    dx = xsize / (xnum - 1)
    dy = ysize / (ynum - 1)
    return min(dx, dy)
end

end # module