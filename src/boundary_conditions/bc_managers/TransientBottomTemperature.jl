module TransientBottomTemperature

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ModelDataContainer: set_parameter!
import EarthBox.InitializationTools: update_use_flags
import EarthBox.ConversionFuncs: celsius_to_kelvin
import EarthBox.ConversionFuncs: kelvin_to_celsius
import EarthBox.ConversionFuncs: seconds_to_years
import EarthBox.ParameterRegistry: get_eb_parameters

import ..SetBoundaryConditions
import ..TransientBottomCalculator

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_bottom_transient::Symbol
    temperature_bottom_transient::Symbol
    start_time_bottom_transient::Symbol
    end_time_bottom_transient::Symbol
    delta_temperature_transient::Symbol
    temperature_base_lith::Symbol
    adiabatic_gradient::Symbol
end

""" 
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize transient bottom temperature boundary conditions by performing
the following calculations:

1. Calculate the distance from the base of the lithosphere to the bottom of the model.
2. Calculate the final cooler temperature at the bottom of the model 
    `temperature_bottom_cooler_final` using the temperature at the base of the 
    lithosphere `temperature_base_lith`, the adiabatic gradient `adiabatic_gradient`, 
    and the distance from the base of the lithosphere to the bottom of the model.
2. Calculate the initial elevated temperature at the base of the lithosphere
    `temperature_base_lith_warmer_initial` using the temperature at the base of the 
    `temperature_base_lith' and the elevated temperature anomaly 
    `delta_temperature_transient`.
3. Calculate the initial elevated temperature at the bottom of the model 
    `temperature_bottom_warmer_initial` using the initial elevated temperature 
    at the base of the lithosphere `temperature_base_lith_warmer_initial`, the 
    distance from the base of the lithosphere to the bottom of the model, and the 
    adiabatic gradient `adiabatic_gradient`.

The transient temperature will be the final cooler temperature at the bottom of the model 
`temperature_bottom_cooler_final` which will be applied between the start and end times 
specified by the user. This is used to approximate the effects of a plume that rapidly
ascends through the lithosphere, spreads out at the base of the lithosphere, and 
cools over time.

# Arguments
- `model::`[`ModelData`](@ref ModelData): Model data object containing model parameters and arrays.

# Optional Keyword Arguments
- `iuse_bottom_transient::Int64`: 
    - $(PDATA.iuse_bottom_transient.description)
- `delta_temperature_transient::Float64`: 
    - $(PDATA.delta_temperature_transient.description)
- `temperature_base_lith::Float64`: 
    - $(PDATA.temperature_base_lith.description)
    - For cases where an transient bottom temperatures this parameter is the final cooler temperature 
       at the base of the lithosphere. If an AnalyticalThreeLayer model is used, the temperature
       at the base of the lithosphere (`temperature_base_lith`) is updated to the warmer value
       when calculating the initial geotherm.
- `start_time_bottom_transient::Float64`: 
    - $(PDATA.start_time_bottom_transient.description)
- `end_time_bottom_transient::Float64`: 
    - $(PDATA.end_time_bottom_transient.description)
- `adiabatic_gradient::Float64`: 
    - $(PDATA.adiabatic_gradient.description)

"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    update_parameters_for_transient_cooling!(model)
    return nothing
end

function update_parameters_for_transient_cooling!(model::ModelData)::Nothing
    temp_params  = model.bcs.parameters.temperature
    geom_params  = model.geometry.parameters
    grids_params = model.grids.parameters
    heat_params  = model.heat_equation.parameters

    adiabatic_gradient          = heat_params.initial_condition.adiabatic_gradient.value
    delta_temperature_transient = temp_params.delta_temperature_transient.value
    temperature_base_lith       = heat_params.steady_state.temperature_base_lith.value
    ysize                       = grids_params.geometry.ysize.value
    thick_air                   = geom_params.sticky_air_geometry.thick_air.value
    thick_lith                  = geom_params.earth_layering.thick_lith.value
    if using_transient_bottom_temperature(model)
        (
            temperature_base_lith_warmer_initial,
            temperature_bottom_warmer_initial, 
            temperature_bottom_cooler_final
        ) = TransientBottomCalculator.calculate_thermal_bcs_for_transient_cooling(
            temperature_base_lith, delta_temperature_transient, adiabatic_gradient,
            ysize, thick_air, thick_lith
            )
        temp_params.temperature_base_lith_warmer_initial.value = temperature_base_lith_warmer_initial
        temp_params.temperature_bottom_warmer_initial.value    = temperature_bottom_warmer_initial
        temp_params.temperature_bottom_cooler_final.value      = temperature_bottom_cooler_final
        temp_params.temperature_bottom.value                   = temperature_bottom_warmer_initial
        temp_params.temperature_bottom_transient.value         = temperature_bottom_cooler_final
        print_info("Calculated transient bottom temperature boundary conditions", level=1)
        print_info("Inputs:", level=2)
        print_info("temperature_base_lith (C): $(kelvin_to_celsius(temperature_base_lith))", level=3)
        print_info("delta_temperature_transient (C): $(delta_temperature_transient)", level=3)
        print_info("adiabatic_gradient (K/km): $(adiabatic_gradient)", level=3)
        print_info("ysize (km): $(ysize)", level=3)
        print_info("thick_air (km): $(thick_air)", level=3)
        print_info("thick_lith (km): $(thick_lith)", level=3)
        print_info("Outputs", level=2)
        print_info("temperature_base_lith_warmer_initial (C): $(kelvin_to_celsius(temperature_base_lith_warmer_initial))", level=3)
        print_info("temperature_bottom_warmer_initial (C): $(kelvin_to_celsius(temperature_bottom_warmer_initial))", level=3)
        print_info("temperature_bottom_cooler_final (C): $(kelvin_to_celsius(temperature_bottom_cooler_final))", level=3)
        print_info("Set Model Data", level=2)
        print_info("temperature_bottom (C): $(kelvin_to_celsius(temp_params.temperature_bottom.value))", level=3)
        print_info("temperature_bottom_transient (C): $(kelvin_to_celsius(temp_params.temperature_bottom_transient.value))", level=3)
    end
    return nothing
end

function using_transient_bottom_temperature(model::ModelData)::Bool
    temp_params  = model.bcs.parameters.temperature
    geom_params  = model.geometry.parameters
    grids_params = model.grids.parameters
    heat_params  = model.heat_equation.parameters

    iuse_bottom_transient       = temp_params.iuse_bottom_transient.value
    adiabatic_gradient          = heat_params.initial_condition.adiabatic_gradient.value
    delta_temperature_transient = temp_params.delta_temperature_transient.value
    temperature_base_lith       = heat_params.steady_state.temperature_base_lith.value
    ysize                       = grids_params.geometry.ysize.value
    thick_air                   = geom_params.sticky_air_geometry.thick_air.value
    thick_lith                  = geom_params.earth_layering.thick_lith.value

    print_info("Parameters for transient bottom temperature model:", level=1)
    print_info("delta_temperature_transient: $(delta_temperature_transient)", level=2)
    print_info("temperature_base_lith: $(temperature_base_lith)", level=2)
    print_info("adiabatic_gradient: $(adiabatic_gradient)", level=2)
    print_info("thick_lith: $(thick_lith)", level=2)
    print_info("thick_air: $(thick_air)", level=2)
    print_info("ysize: $(ysize)", level=2)
    print_info("iuse_bottom_transient: $(iuse_bottom_transient)", level=2)
    print_info("", level=1)

    iuse_bottom_transient = model.bcs.parameters.temperature.iuse_bottom_transient.value
    if (
        iuse_bottom_transient == 1 && !isnan(delta_temperature_transient)
        && temperature_base_lith > 0.0 && adiabatic_gradient > 0.0 
        && thick_lith > 0.0 && thick_air > 0.0 && ysize > 0.0
    )
        return true
    end
    return false
end

function update_bottom_temperature!(model::ModelData)::Nothing
    iuse_bottom_transient = model.bcs.parameters.temperature.iuse_bottom_transient.value
    if iuse_bottom_transient == 1
        apply_transient_bottom_temperature!(model)
    end
end

function apply_transient_bottom_temperature!(model::ModelData)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_myr = seconds_to_years(timesum)/1e6
    start_myr = model.bcs.parameters.temperature.start_time_bottom_transient.value
    end_myr = model.bcs.parameters.temperature.end_time_bottom_transient.value

    temperature = model.bcs.parameters.temperature
    temperature_bottom_original = temperature.temperature_bottom_original.value
    temperature_bottom = temperature.temperature_bottom.value

    if temperature_bottom_original < 0.0
        set_original_bottom_temperature!(model, temperature_bottom)
    end

    if start_myr <= timesum_myr <= end_myr
        temperature_bottom = model.bcs.parameters.temperature.temperature_bottom_transient.value
        print_info("Bottom temperature (Transient) (C): $(kelvin_to_celsius(temperature_bottom))", level=2)
    else
        temperature_bottom = model.bcs.parameters.temperature.temperature_bottom_original.value
        print_info("Bottom temperature (Original) (C): $(kelvin_to_celsius(temperature_bottom))", level=2)
    end

    set_bottom_temperature!(model, temperature_bottom)
    SetBoundaryConditions.BottomBoundary.temperature!(model)
end

function set_original_bottom_temperature!(model::ModelData, temperature_bottom::Float64)
    temperature = model.bcs.parameters.temperature
    temperature.temperature_bottom_original.value = temperature_bottom
end

function set_bottom_temperature!(model::ModelData, temperature_bottom_new::Float64)
    temperature = model.bcs.parameters.temperature
    temperature.temperature_bottom.value = temperature_bottom_new
end

function using_bottom_transient(model::ModelData)::Bool
    return model.bcs.parameters.temperature.iuse_bottom_transient.value == 1
end

end # module 