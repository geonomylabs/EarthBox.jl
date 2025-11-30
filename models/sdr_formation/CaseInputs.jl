"""
    CaseInputs

Custom module for defining case inputs for SDR formation model.

To print a table of case inputs from he Julia REPL, run:

```julia
using EarthBox
PRINT_SETTINGS.print_case_info = true;
include("CaseInputs.jl");
CaseInputs.define_case_inputs();
```
"""
module CaseInputs

using EarthBox

const PARAMS = get_eb_parameters()

function define_case_parameters(model_case_name::String="case0")::CaseType
    base_case = CaseType(
        # ThermalBottomTransientBoundaryCondition
        PARAMS.delta_temperature_transient.name => CaseParameter(100.0, "deltaK"),
        # GlobalPlasticityLoop
        PARAMS.tolerance_picard.name => CaseParameter(1e-3, "None"),
        PARAMS.nglobal.name          => CaseParameter(100, "None"),
        # AdvectionModel
        PARAMS.marker_cell_displ_max.name => CaseParameter(0.3, "fraction"),
        # MeltDamageModel
        PARAMS.melt_damage_distance.name => CaseParameter(2_500.0, "m"),
        PARAMS.melt_damage_factor.name   => CaseParameter(10.0, "None"),
        # Extrusion Model
        PARAMS.characteristic_flow_length_subaerial.name => CaseParameter(20_000.0, "m"),
        PARAMS.eruption_interval_yr.name                 => CaseParameter(50_000.0, "yr"),
        PARAMS.extrusion_volume_factor_max.name           => CaseParameter(0.5, "None"),
        # Sediment Transport Model
        PARAMS.subaerial_transport_coefficient.name => CaseParameter(1.0e-4, "None"),
        PARAMS.pelagic_sedimentation_rate.name      => CaseParameter(0.0, "mm/yr"),
        # Material Overrides (see [Material Override Using CaseType](@ref))
        PARAMS.latent_heat_mantle.name        => CaseParameter(400_000.0, "J/kg"),
        PARAMS.latent_heat_oceanic_crust.name => CaseParameter(400_000.0, "J/kg"),
        PARAMS.mantle_solidus.name            => CaseParameter(solidus_names.PeridotiteKatz2003, "None"),
        PARAMS.mantle_liquidus.name           => CaseParameter(liquidus_names.PeridotiteKatz2003, "None"),
    )
    # Initialize using the base case
    case_inputs = initialize_cases(base_case)
    # Update case_inputs with 4 cases with variable delta temperature transient
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = 1,
        parameter_name = PARAMS.delta_temperature_transient.name,
        values         = [50.0, 75.0, 125.0, 150.0],
        units          = "deltaK"
    )
    # Update case_inputs with 4 cases with variable characteristic flow length subaerial
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.characteristic_flow_length_submarine.name,
        values         = [10_000.0, 40_000.0, 80_000.0, 100_000.0],
        units          = "m"
    )
    # Update case_inputs with 4 cases with variable melt damage factor
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.melt_damage_factor.name,
        values         = [1.0, 1.25, 2.5, 5.0],
        units          = "None"
    )
    # Update case_inputs with 3 cases with variable melt damage distance
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.melt_damage_distance.name,
        values         = [625.0, 1_250.0, 5_000.0],
        units          = "m"
    )
    # Update case_inputs with 3 cases with variable eruption duration
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.eruption_interval_yr.name,
        values         = [12_500.0, 25_000.0, 100_000.0],
        units          = "yr"
    )
    # Update case_inputs with 6 cases with variable maximum extrusion volume factor @ eruption duration = 50_000 yrs
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.extrusion_volume_factor_max.name,
        values               = [0.1, 0.2, 0.3, 0.4, 0.6, 0.7],
        units                = "None",
        fixed_parameter_name = PARAMS.eruption_interval_yr.name,
        fixed_value          = 50_000.0,
        fixed_units          = "yr"
    )
    # Update case_inputs with 4 cases with variable maximum extrusion volume factor @ eruption duration = 75_000 yrs
    case_id = define_case_group!(
        case_inputs,
        case_id_ini            = case_id + 1,
        parameter_name         = PARAMS.extrusion_volume_factor_max.name,
        values                 = [0.1, 0.3, 0.5, 0.7],
        units                  = "None",
        fixed_parameter_name   = PARAMS.eruption_interval_yr.name,
        fixed_value            = 75_000.0,
        fixed_units            = "yr"
    )
    # Update case_inputs with 4 cases with variable maximum extrusion volume factor @ eruption duration = 100_000 yrs
    case_id = define_case_group!(
        case_inputs,
        case_id_ini            = case_id + 1,
        parameter_name         = PARAMS.extrusion_volume_factor_max.name,
        values                 = [0.1, 0.2, 0.3, 0.7],
        units                  = "None",
        fixed_parameter_name   = PARAMS.eruption_interval_yr.name,
        fixed_value            = 100_000.0,
        fixed_units            = "yr"
    )
    # Update case_inputs with 4 cases with variable maximum extrusion volume factor @ damage factor = 1
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.extrusion_volume_factor_max.name,
        values               = [0.1, 0.2, 0.3, 0.4],
        units                = "None",
        fixed_parameter_name = PARAMS.melt_damage_factor.name,
        fixed_value          = 1.0,
        fixed_units          = "None"
    )
    print_case_info(
        case_inputs=case_inputs, 
        case_id_max=case_id, 
        target_names=[
            PARAMS.delta_temperature_transient.name,
            PARAMS.extrusion_volume_factor.name,
            PARAMS.extrusion_volume_factor_max.name,
            PARAMS.iuse_melt_damage.name,
            PARAMS.melt_damage_factor.name,
            PARAMS.melt_damage_distance.name,
            PARAMS.characteristic_flow_length_subaerial.name,
            PARAMS.eruption_interval_yr.name,
            PARAMS.tolerance_picard.name,
            PARAMS.nglobal.name,
            PARAMS.marker_cell_displ_max.name,
            PARAMS.subaerial_transport_coefficient.name,
            PARAMS.pelagic_sedimentation_rate.name
        ]
    )
    case_parameters = case_inputs[model_case_name]
    convert_case_parameters_to_standard_units!(case_parameters)
    return case_parameters
end

end # module
