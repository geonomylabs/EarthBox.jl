module HeatSolver

include("rhs/HeatRhs.jl")
include("solver/SystemSolver.jl")
include("heat_loop/HeatLoop.jl")
include("subgrid_heat/SubgridDiffusion.jl")
include("marker_temperature/MarkerTemperature.jl")

using InteractiveUtils
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState
import EarthBox.Markers.MarkerTemperature.InitManager: PlumeUpdate
import EarthBox.Arrays: ArrayUtils
import .HeatLoop
import .HeatRhs
import .SubgridDiffusion
import .MarkerTemperature

const PDATA = get_eb_parameters()

struct ValidInputNames
    max_temp_change::Symbol
    iuse_heat::Symbol
    iuse_adiabatic_heating::Symbol
    iuse_shear_heating::Symbol
    iuse_sticky_correction::Symbol
end

"""
    initialize!(
        model::ModelData;
        kwargs...
    )::Nothing

Initialize heat solver parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing 
    model parameters and arrays.

# Keyword Arguments
- `$(PDATA.max_temp_change.name)::Float64`:
    - $(PDATA.max_temp_change.description)
- `$(PDATA.iuse_heat.name)::Int64`:
    - $(PDATA.iuse_heat.description)
- `$(PDATA.iuse_adiabatic_heating.name)::Int64`:
    - $(PDATA.iuse_adiabatic_heating.description)
- `$(PDATA.iuse_shear_heating.name)::Int64`:
    - $(PDATA.iuse_shear_heating.description)
- `$(PDATA.iuse_sticky_correction.name)::Int64`:
    - $(PDATA.iuse_sticky_correction.description)

"""
function initialize!(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Heat Solver Option")
end

function get_option_id_from_model(model::ModelData)::Int
    return model.heat_equation.parameters.heat_options.itype_heat.value
end

function get_stype_from_model(model::ModelData)::String
    return model.heat_equation.parameters.heat_options.stype_heat.value
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.heat_equation.parameters.heat_options.itype_heat.value = option_id
    return nothing
end

""" Update marker temperature for grid and subgrid changes.

This function orchestrates steps that update marker temperature
taking into account subgrid temperature changes and conductive 
temperature changes on basic grid nodes.

# Arguments
- `model::ModelData`: The model data.
- `solver_config::SolverConfig`: The solver configuration.

"""
function update_marker_temp_for_grid_and_subgrid_changes!(
    model::ModelData,
    solver_config::SolverConfigState,
    inside_flags::Vector{Int8}
)::Nothing
    if use_heat_solver(model)
        calculate_heat_rhs_terms_grid!(model)
        initialize_old_temperature_array!(model)
        @timeit_memit "Finished conductive grid heat solver loop" begin
            HeatLoop.run_conductive_heat_solver_loop!(model, solver_config)
        end
        @timeit_memit "Finished calculating temperature change on basic grid" begin
            calculate_grid_temp_change!(model)
        end
        if use_subgrid_diffusion(model)
            @timeit_memit "Finished calculating subgrid temperature changes" begin
                SubgridDiffusion.subgrid_heat_diffusion_update!(model, inside_flags)
            end
            @timeit_memit "Finished removing subgrid thermal effects from grid changes" begin
                SubgridDiffusion.remove_subgrid_effects_from_grid_changes!(model)
            end
        end
        @timeit_memit "Finished adding remaining grid change to marker temperature" begin
            MarkerTemperature.add_remaining_grid_change_to_marker_temperature!(
                model, inside_flags)
        end
    end
    return nothing
end

""" Calculate right-hand-side terms for heat equation on basic grid.

This method calls functions that calculate right-hand terms array
for the heat equation on the basic grid. Radiogenic heat 
production is included using the radiogenic heat production 
array. Adiabatic heating terms are calculated using the adiabatic 
coefficients arrays, which is the product of expansivity and 
temperature, the density array, and the staggered velocity 
solution arrays. Shear heating is included using visco-elasto-plastic 
stress arrays and viscoplastic viscosity arrays.

# Arguments
- `model::ModelData`: The model data.

"""
function calculate_heat_rhs_terms_grid!(model::ModelData)::Nothing
    @timeit_memit "Finished calculating heat right-hand-side terms" begin
        HeatRhs.initialize_rhs_terms_array_for_heat!(model)
        if model.heat_equation.parameters.heat_options.iuse_adiabatic_heating.value == 1
            HeatRhs.add_adiabatic_terms_to_rhs_grid!(model)
        end
        if model.heat_equation.parameters.heat_options.iuse_shear_heating.value == 1
            HeatRhs.add_friction_terms_to_rhs_grid!(model)
        end
    end
    return nothing
end

function use_heat_solver(model::ModelData)::Bool
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    max_temp_change = model.heat_equation.parameters.temp_change_limit.max_temp_change.value
    iuse_heat = model.heat_equation.parameters.heat_options.iuse_heat.value

    if timestep > 0 && max_temp_change > 0 && iuse_heat == 1
        return true
    end
    return false
end

function use_subgrid_diffusion(model::ModelData)::Bool
    (
        subgrid_diff_coef_temp
    ) = model.markers.parameters.subgrid_diffusion.subgrid_diff_coef_temp.value
    if subgrid_diff_coef_temp > 0
        return true
    end
    return false
end

function apply_sticky_temperature_correction!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    if model.heat_equation.parameters.heat_options.iuse_sticky_correction.value == 1
        @timeit_memit "Finished applying sticky temperature correction" begin
            MarkerTemperature.apply_sticky_temperature_correction!(model, inside_flags)
        end
    end
    return nothing
end

function manage_plume_injection!(model::ModelData)::Nothing
    if model.geometry.parameters.plume.iuse_plume.value == 1
        @timeit_memit "Finished managing plume injection" begin
            PlumeUpdate.inject_plume!(model)
        end
    end
    return nothing
end

""" Initialize the old temperature solution array.
This function initializes the old temperature solution array by
copying the basic grid temperature array that was previously
interpolated from advected markers.
"""
function initialize_old_temperature_array!(model::ModelData)::Nothing
    model.heat_equation.arrays.temperature.tk0.array .= 
        model.heat_equation.arrays.temperature.tk1.array
    return nothing
end

""" Calculate grid temperature change.
This function calculates the temperature change array on the basic 
grid by subtracting the temperature transport array (interpolated 
from advected markers) from the new temperature solution array based
on the current heat conduction solution.
"""
function calculate_grid_temp_change!(model::ModelData)::Nothing
    model.heat_equation.arrays.temperature.dtk1.array .= 
        model.heat_equation.arrays.temperature.tk2.array .- 
        model.heat_equation.arrays.temperature.tk1.array
    return nothing
end

end # module 