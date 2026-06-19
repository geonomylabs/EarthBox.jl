"""
    CaseInputs

Custom module for defining case inputs.

| Case   | delta_     | extrusion_ | melt_      | maximum_   | dike_      | characteri | eruption_  | subaerial_ | full_      |
|        | temperatur | volume_    | damage_    | damage_    | fluid_     | stic_flow_ | interval_  | transport_ | velocity_  |
|        | e_transien | factor_max | factor     | probabilit | marker_    | length_    | yr         | coefficien | extension  |
|        | t          |            |            | y          | fraction   | subaerial  |            | t          |            |
|--------|------------|------------|------------|------------|------------|------------|------------|------------|------------|
| 0      | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 1      | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 5.00e+04   | 1.00e-04   | 2          |
| 2      | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 1.00e+05   | 1.00e-04   | 2          |
| 3      | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 4.00e+05   | 1.00e-04   | 2          |
| 4      | 25         | 0.5        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 5      | 50         | 0.5        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 6      | 150        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 7      | 100        | 0.5        | 10         | 1          | 0          | 2.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 8      | 100        | 0.5        | 10         | 1          | 0          | 6.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 9      | 100        | 0.5        | 10         | 1          | 0          | 8.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 10     | 100        | 0.3        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 11     | 100        | 0.4        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 12     | 100        | 0.6        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 13     | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 0.01       | 2          |
| 14     | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 1          |
| 15     | 100        | 0.5        | 10         | 1          | 1          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 16     | 100        | 0.5        | 10         | 0.75       | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 17     | 100        | 0.5        | 1          | 1          | 0          | 4.00e+04   | 5.00e+04   | 1.00e-04   | 2          |
| 18     | 100        | 0.5        | 1          | 1          | 0          | 4.00e+04   | 1.00e+05   | 1.00e-04   | 2          |
| 19     | 100        | 0.5        | 1          | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 20     | 100        | 0.5        | 1          | 1          | 0          | 2.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 21     | 100        | 0.5        | 1          | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 22     | 100        | 0.5        | 1          | 1          | 0          | 6.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 23     | 25         | 0.5        | 1          | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 24     | 50         | 0.5        | 1          | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 25     | 150        | 0.5        | 1          | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 26     | 100        | 0.3        | 1          | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 27     | 100        | 0.4        | 1          | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 28     | 100        | 0.6        | 1          | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 29     | 100        | 0.5        | 1          | 1          | 0          | 8.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 30     | 100        | 0.5        | 1          | 1          | 0          | 4.00e+04   | 4.00e+05   | 1.00e-04   | 2          |
| 31     | 150        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 1          |
| 32     | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 4.00e+05   | 1.00e-04   | 1          |
| 33     | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 0.01       | 1          |
| 34     | 150        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 4.00e+05   | 1.00e-04   | 1          |
| 35     | 150        | 0.5        | 10         | 1          | 0          | 6.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 36     | 150        | 0.5        | 10         | 1          | 0          | 8.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 37     | 150        | 0.3        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 38     | 150        | 0.6        | 10         | 1          | 0          | 4.00e+04   | 2.00e+05   | 1.00e-04   | 2          |
| 39     | 150        | 0.5        | 10         | 1          | 0          | 8.00e+04   | 4.00e+05   | 1.00e-04   | 1          |
| 40     | 150        | 0.5        | 10         | 1          | 0          | 8.00e+04   | 4.00e+05   | 0.01       | 1          |
| 41     | 100        | 0.5        | 10         | 1          | 0          | 4.00e+04   | 4.00e+05   | 0.01       | 1          |
| 42     | 100        | 0.5        | 10         | 1          | 0          | 8.00e+04   | 4.00e+05   | 1.00e-04   | 1          |

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
        # Velocity boundary condition
        PARAMS.full_velocity_extension.name => CaseParameter(2.0, "cm/yr"),
        # ThermalBottomTransientBoundaryCondition
        PARAMS.delta_temperature_transient.name => CaseParameter(100.0, "deltaK"),
        # MeltDamageModel
        PARAMS.melt_damage_distance.name => CaseParameter(2_500.0, "m"),
        PARAMS.melt_damage_factor.name   => CaseParameter(10.0, "None"),
        PARAMS.maximum_damage_probability.name => CaseParameter(1.0, "fraction"),
        PARAMS.dike_fluid_marker_fraction.name => CaseParameter(0.0, "fraction"),
        PARAMS.density_dike_fluid.name => CaseParameter(2_750.0, "kg/m^3"),
        # Extrusion Model
        PARAMS.characteristic_flow_length_subaerial.name => CaseParameter(40_000.0, "m"),
        PARAMS.eruption_interval_yr.name                 => CaseParameter(200_000.0, "yr"),
        PARAMS.extrusion_volume_factor_max.name          => CaseParameter(0.5, "None"),
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
    #
    #**************************** Melt Damage On ****************************
    #
    #****************************
    # Variable eruption intervals
    #****************************
    # case1, case2, case3
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = 1,
        parameter_name = PARAMS.eruption_interval_yr.name,
        values         = [50_000.0, 100_000.0, 400_000.0],
        units          = "yr",
    )
    #*****************
    # Variable delta T
    #*****************
    # case4, case5, case6
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.delta_temperature_transient.name,
        values         = [25.0, 50.0, 150.0],
        units          = "deltaK",
    )
    #*********************
    # Variable flow length
    #*********************
    # case7, case8, case9
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.characteristic_flow_length_subaerial.name,
        values         = [20_000.0, 60_000.0, 80_000.0],
        units          = "m",
    )
    #*********************************
    # Variable extrusion volume factor
    #*********************************
    # case10, case11, case12
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.extrusion_volume_factor_max.name,
        values         = [0.3, 0.4, 0.6],
        units          = "None",
    )
    #*****************************************
    # High subaerial_transport_coefficient
    #*****************************************
    # case13
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.subaerial_transport_coefficient.name,
        values         = [1.0e-2],
        units          = "None",
    )
    #*****************************************
    # Low full_velocity_extension
    #*****************************************
    # case14
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
    )
    #******************************
    # Dike fluid marker fraction on
    #******************************
    # case15
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.dike_fluid_marker_fraction.name,
        values         = [1.0],
        units          = "fraction",
    )
    #******************************
    # Low maximum damage probability
    #******************************
    # case16
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.maximum_damage_probability.name,
        values         = [0.75],
        units          = "fraction",
    )
    #
    #**************************** Melt Damage Off ****************************
    #
    #****************************
    # Variable eruption intervals
    #****************************
    # case17, case18, case19
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.eruption_interval_yr.name,
        values         = [50_000.0, 100_000.0, 200_000.0],
        units          = "yr",
        fixed                = [
            (PARAMS.melt_damage_factor.name, 1.0, "None"),
            ]
    )
    #*********************
    # Variable flow length
    #*********************
    # case20, case21, case22
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.characteristic_flow_length_subaerial.name,
        values         = [20_000.0, 40_000.0, 60_000.0],
        units          = "m",
        fixed                = [
            (PARAMS.melt_damage_factor.name, 1.0, "None"),
            ]
    )
    #*************************
    # Variable delta T
    #*************************
    # case23, case24, case25
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.delta_temperature_transient.name,
        values         = [25.0, 50.0, 150.0],
        units          = "deltaK",
        fixed                = [
            (PARAMS.melt_damage_factor.name, 1.0, "None"),
            ]
    )
    #*************************
    # Variable extrusion volume factor
    #*************************
    # case26, case27, case28
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.extrusion_volume_factor_max.name,
        values         = [0.3, 0.4, 0.6],
        units          = "None",
        fixed                = [
            (PARAMS.melt_damage_factor.name, 1.0, "None"),
            ]
    )
    #******************
    # Additional cases 
    #******************
    # case29
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.characteristic_flow_length_subaerial.name,
        values         = [80_000.0],
        units          = "m",
        fixed                = [
            (PARAMS.melt_damage_factor.name, 1.0, "None"),
            ]
    )
    # case30
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.eruption_interval_yr.name,
        values         = [400_000.0],
        units          = "yr",
        fixed                = [
            (PARAMS.melt_damage_factor.name, 1.0, "None"),
            ]
    )
    #
    #****** Additional cases with low full_velocity_extension and melt damage on *******
    #
    # case31
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
        fixed                = [
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            ]
    )
    # case32
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
        fixed                = [
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            ]
    )
    # case33
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
        fixed                = [
            (PARAMS.subaerial_transport_coefficient.name, 1.0e-2, "None"),
            ]
    )
    # case34
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
        fixed                = [
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            ]
    )
    #
    #****** Additional cases with delta T = 150C *******
    #
    #*********************
    # Variable flow length
    #*********************
    # case35, case36
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.characteristic_flow_length_subaerial.name,
        values         = [60_000.0, 80_000.0],
        units          = "m",
        fixed                = [
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            ]
    )
    #*********************************
    # Variable extrusion volume factor
    #*********************************
    # case37, case38
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.extrusion_volume_factor_max.name,
        values         = [0.3, 0.6],
        units          = "None",
        fixed                = [
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            ]
    )
    #
    #************* Additional cases with velocity = 1.0 cm/yr *************
    #
    # case39 (similar to case 34 but with longer flow length)
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
        fixed                = [
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            (PARAMS.characteristic_flow_length_subaerial.name, 80_000.0, "m"),
            ]
    )
    # case40 (similar to case 34 but with longer flow length and greater subaerial transport)
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
        fixed                = [
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            (PARAMS.delta_temperature_transient.name, 150.0, "deltaK"),
            (PARAMS.characteristic_flow_length_subaerial.name, 80_000.0, "m"),
            (PARAMS.subaerial_transport_coefficient.name, 1.0e-2, "None"),
            ]
    )
    # case 41 (similar to case 32 but with greater subaerial transport)
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
        fixed                = [
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            (PARAMS.subaerial_transport_coefficient.name, 1.0e-2, "None"),
            ]
    )
    # case 42 (similar to case 32 but with longer flow length)
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.full_velocity_extension.name,
        values         = [1.0],
        units          = "cm/yr",
        fixed                = [
            (PARAMS.eruption_interval_yr.name, 400_000.0, "yr"),
            (PARAMS.characteristic_flow_length_subaerial.name, 80_000.0, "m"),
            ]
    )

    print_case_info(
        case_inputs=case_inputs, 
        case_id_max=case_id, 
        target_names=[
            PARAMS.delta_temperature_transient.name,
            PARAMS.extrusion_volume_factor_max.name,
            PARAMS.melt_damage_factor.name,
            PARAMS.maximum_damage_probability.name,
            PARAMS.dike_fluid_marker_fraction.name,
            PARAMS.characteristic_flow_length_subaerial.name,
            PARAMS.eruption_interval_yr.name,
            PARAMS.subaerial_transport_coefficient.name,
            PARAMS.full_velocity_extension.name,
        ]
    )
    case_parameters = case_inputs[model_case_name]
    convert_case_parameters_to_standard_units!(case_parameters)
    return case_parameters
end

end # module
