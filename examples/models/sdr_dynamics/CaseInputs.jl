"""
    CaseInputs

Custom module for defining case inputs.

| Case   | delta_     | extrusion_ | melt_      | melt_      | maximum_   | dike_      | characteri | eruption_  | Node
|        | temperatur | volume_    | damage_    | damage_    | damage_    | fluid_     | stic_flow_ | interval_  |
|        | e_transien | factor_max | factor     | distance   | probabilit | marker_    | length_    | yr         |
|        | t          |            |            |            | y          | fraction   | subaerial  |            | 
|--------|------------|------------|------------|------------|------------|------------|------------|------------|------------|
| 0      | 100        | 0.5        | 10         | 2500       | 1          | 1          | 4.00e+04   | 2.00e+05   | peridotite (running)
| 1      | 100        | 0.5        | 1          | 2500       | 1          | 0          | 4.00e+04   | 2.00e+05   | peridotite (running)
| 2      | 100        | 0.5        | 10         | 2500       | 1          | 1          | 6.00e+04   | 2.00e+05   | eclogite (running)
| 3      | 100        | 0.5        | 10         | 2500       | 1          | 1          | 4.00e+04   | 4.00e+05   | eclogite (running)
| 4      | 100        | 0.5        | 10         | 2500       | 1          | 1          | 6.00e+04   | 4.00e+05   | lherzolite (running)
| 5      | 100        | 0.7        | 10         | 2500       | 1          | 1          | 4.00e+04   | 2.00e+05   | lherzolite (running)
| 6      | 100        | 0.7        | 10         | 2500       | 1          | 1          | 6.00e+04   | 2.00e+05   | plagioclase (running)
| 7      | 100        | 0.7        | 10         | 2500       | 1          | 1          | 4.00e+04   | 4.00e+05   | plagioclase (running)
| 8      | 100        | 0.7        | 10         | 2500       | 1          | 1          | 6.00e+04   | 4.00e+05   | chromite (running)
| 9      | 150        | 0.5        | 10         | 2500       | 1          | 1          | 4.00e+04   | 2.00e+05   | chromite (running)
| 10     | 150        | 0.5        | 10         | 2500       | 1          | 1          | 6.00e+04   | 2.00e+05   | pyroxene (running)
| 11     | 150        | 0.5        | 10         | 2500       | 1          | 1          | 4.00e+04   | 4.00e+05   | pyroxene (running)
| 12     | 150        | 0.5        | 10         | 2500       | 1          | 1          | 6.00e+04   | 4.00e+05   | spinel (running)
| 13     | 150        | 0.7        | 10         | 2500       | 1          | 1          | 4.00e+04   | 2.00e+05   | spinel (running)
| 14     | 150        | 0.7        | 10         | 2500       | 1          | 1          | 6.00e+04   | 2.00e+05   | komatiite (running)
| 15     | 150        | 0.7        | 10         | 2500       | 1          | 1          | 4.00e+04   | 4.00e+05   | komatiite (running)
| 16     | 150        | 0.7        | 10         | 2500       | 1          | 1          | 6.00e+04   | 4.00e+05   | dunite (running)
| 17     | 100        | 0.5        | 10         | 2500       | 1          | 1          | 4.00e+04   | 5.00e+04   |
| 18     | 100        | 0.5        | 10         | 2500       | 1          | 1          | 4.00e+04   | 1.00e+05   |
| 19     | 100        | 0.5        | 1          | 2500       | 1          | 0          | 4.00e+04   | 5.00e+04   |
| 20     | 100        | 0.5        | 1          | 2500       | 1          | 0          | 4.00e+04   | 1.00e+05   |
| 21     | 100        | 0.5        | 10         | 2500       | 1          | 0          | 4.00e+04   | 5.00e+04   |
| 22     | 100        | 0.5        | 10         | 2500       | 1          | 0          | 4.00e+04   | 1.00e+05   |
| 23     | 100        | 0.5        | 10         | 2500       | 1          | 0          | 4.00e+04   | 2.00e+05   | spinel (running)

To print a table of case inputs from he Julia REPL, run:

```julia
using EarthBox
PRINT_SETTINGS.print_case_info = true;
include("CaseInputs.jl");
CaseInputs.define_case_parameters();
```
"""
module CaseInputs

using EarthBox

const PARAMS = get_eb_parameters()

function define_case_parameters(model_case_name::String="case0")::CaseType
    base_case = CaseType(
        # ThermalBottomTransientBoundaryCondition
        PARAMS.delta_temperature_transient.name => CaseParameter(100.0, "deltaK"),
        # MeltDamageModel
        PARAMS.melt_damage_distance.name => CaseParameter(2_500.0, "m"),
        PARAMS.melt_damage_factor.name   => CaseParameter(10.0, "None"),
        PARAMS.maximum_damage_probability.name => CaseParameter(1.0, "fraction"),
        PARAMS.dike_fluid_marker_fraction.name => CaseParameter(1.0, "fraction"),
        PARAMS.density_dike_fluid.name => CaseParameter(2_750.0, "kg/m^3"),
        # Extrusion Model
        PARAMS.characteristic_flow_length_subaerial.name => CaseParameter(40_000.0, "m"),
        PARAMS.eruption_interval_yr.name                 => CaseParameter(200_000.0, "yr"),
        PARAMS.extrusion_volume_factor_max.name           => CaseParameter(0.5, "None"),
        # Sediment Transport Model
        PARAMS.subaerial_transport_coefficient.name => CaseParameter(1.0e-4, "None"),
        # Material Overrides (see [Material Override Using CaseType](@ref))
        PARAMS.latent_heat_mantle.name        => CaseParameter(400_000.0, "J/kg"),
        PARAMS.latent_heat_oceanic_crust.name => CaseParameter(400_000.0, "J/kg"),
        PARAMS.mantle_solidus.name            => CaseParameter(solidus_names.PeridotiteKatz2003, "None"),
        PARAMS.mantle_liquidus.name           => CaseParameter(liquidus_names.PeridotiteKatz2003, "None"),
        # Rheology
        PARAMS.friction_angle_initial_solidified_basalt.name => CaseParameter(30.0, "degrees"),
        PARAMS.friction_angle_final_solidified_basalt.name => CaseParameter(7.0, "degrees"),
    )
    # Initialize using the base case
    case_inputs = initialize_cases(base_case)
    #***************
    # No melt damage
    #***************
    # Case without melt damage and no density reduction from dike fluid
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = 1,
        parameter_name = PARAMS.melt_damage_factor.name,
        values         = [1.0],
        units          = "None",
        fixed                = [
            (PARAMS.dike_fluid_marker_fraction.name, 0.0, "fraction"),
            ]
    )
    #**********************************
    # Delta T = 100C, extrusion at 50 %
    #**********************************
    # Base case but with longer flow length
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.characteristic_flow_length_subaerial.name,
        values         = [60_000.0],
        units          = "m"
    )
    # 400kyr eruption interval
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.characteristic_flow_length_subaerial.name,
        values               = [40_000.0, 60_000.0],
        units                = "m",
        fixed                = [
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            ]
    )
    #**********************************
    # Delta T = 100C, extrusion at 70 %
    #**********************************
    # Eruption duration of 200 Kyr
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.characteristic_flow_length_subaerial.name,
        values               = [40_000.0, 60_000.0],
        units                = "m",
        fixed                = [
            (PARAMS.extrusion_volume_factor_max.name, 0.7, "None"),
            ]
    )
    # Eruption duration of 400 Kyr
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.characteristic_flow_length_subaerial.name,
        values               = [40_000.0, 60_000.0],
        units                = "m",
        fixed                = [
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            (PARAMS.extrusion_volume_factor_max.name, 0.7, "None"),
            ]
    )
    #**********************************
    # Delta T = 150C, extrusion at 50 %
    #**********************************
    # Eruption duration of 200 Kyr
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.characteristic_flow_length_subaerial.name,
        values               = [40_000.0, 60_000.0],
        units                = "m",
        fixed                = [
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            ]
    )
    # Eruption duration of 400 Kyr
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.characteristic_flow_length_subaerial.name,
        values               = [40_000.0, 60_000.0],
        units                = "m",
        fixed                = [
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            ]
    )
    #**********************************
    # Delta T = 150C, extrusion at 70 %
    #**********************************
    # Eruption duration of 200 Kyr
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.characteristic_flow_length_subaerial.name,
        values               = [40_000.0, 60_000.0],
        units                = "m",
        fixed                = [
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            (PARAMS.extrusion_volume_factor_max.name, 0.7, "None"),
            ]
    )
    # Eruption duration of 400 Kyr
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.characteristic_flow_length_subaerial.name,
        values               = [40_000.0, 60_000.0],
        units                = "m",
        fixed                = [
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            (PARAMS.extrusion_volume_factor_max.name, 0.7, "None"),
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            ]
    )
    #***************************************************************
    # Delta T = 100C, extrusion at 50 %, shorter eruption intervals
    #***************************************************************
    # Base case but with longer flow length
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.eruption_interval_yr.name,
        values         = [50_000.0, 100_000.0],
        units          = "yr"
    )
    #*******************************************
    # No melt damage, shorter eruption intervals
    #*******************************************
    # Case without melt damage and no density reduction from dike fluid
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.eruption_interval_yr.name,
        values         = [50_000.0, 100_000.0],
        units          = "yr",
        fixed                = [
            (PARAMS.melt_damage_factor.name, 1.0, "None"),
            (PARAMS.dike_fluid_marker_fraction.name, 0.0, "fraction"),
            ]
    )
    #*******************************************************
    # No dike density reduction, variable eruption intervals
    #*******************************************************
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.eruption_interval_yr.name,
        values         = [50_000.0, 100_000.0, 200_000.0],
        units          = "yr",
        fixed                = [
            (PARAMS.dike_fluid_marker_fraction.name, 0.0, "fraction"),
            ]
    )
    print_case_info(
        case_inputs=case_inputs, 
        case_id_max=case_id, 
        target_names=[
            PARAMS.delta_temperature_transient.name,
            PARAMS.extrusion_volume_factor_max.name,
            PARAMS.melt_damage_factor.name,
            PARAMS.melt_damage_distance.name,
            PARAMS.maximum_damage_probability.name,
            PARAMS.dike_fluid_marker_fraction.name,
            PARAMS.characteristic_flow_length_subaerial.name,
            PARAMS.eruption_interval_yr.name,
        ]
    )
    case_parameters = case_inputs[model_case_name]
    convert_case_parameters_to_standard_units!(case_parameters)
    return case_parameters
end

end # module
