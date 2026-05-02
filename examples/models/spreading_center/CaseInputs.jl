"""
    CaseInputs

Custom module for defining case inputs for SDR formation model.

| Case   | delta_     | characteri | eruption_  | extrusion_ | maximum_   | melt_      | melt_      | dike_      | density_   | iuse_melt_ |
|        | temperatur | stic_flow_ | interval_  | volume_    | damage_    | damage_    | damage_    | fluid_     | dike_fluid | compaction |
|        | e_transien | length_    | yr         | factor_max | probabilit | distance   | factor     | marker_    |            |            |
|        | t          | subaerial  |            |            | y          |            |            | fraction   |            |            |
|--------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|
| 0      | 0          | 2.00e+04   | 2.00e+05   | 0.5        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 1      | 0          | 2.00e+04   | 1.00e+05   | 0.5        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 2      | 0          | 2.00e+04   | 3.00e+05   | 0.5        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 3      | 0          | 3.00e+04   | 1.00e+05   | 0.5        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 4      | 0          | 3.00e+04   | 3.00e+05   | 0.5        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 5      | 0          | 2.00e+04   | 2.00e+05   | 0.7        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 6      | 0          | 2.00e+04   | 3.00e+05   | 0.7        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 7      | 100        | 2.00e+04   | 2.00e+05   | 0.5        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 8      | 100        | 3.00e+04   | 2.00e+05   | 0.5        | 1          | 2500       | 10         | 1          | 2750       | 0          |
| 9      | 100        | 4.00e+04   | 2.00e+05   | 0.5        | 1          | 2500       | 10         | 1          | 2750       | 0          |

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
        PARAMS.delta_temperature_transient.name => CaseParameter(0.0, "deltaK"),
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
        # Rheology
        PARAMS.friction_angle_initial_solidified_basalt.name => CaseParameter(10.0, "degrees"),
        PARAMS.friction_angle_final_solidified_basalt.name => CaseParameter(1.0, "degrees"),
    )
    # Initialize using the base case
    case_inputs = initialize_cases(base_case)
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = 1,
        parameter_name       = PARAMS.eruption_interval_yr.name,
        values               = [100_000.0, 300_000.0],
        units                = "yr",
    )
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.eruption_interval_yr.name,
        values               = [100_000.0, 300_000.0],
        units                = "yr",
        fixed                = [
            (PARAMS.characteristic_flow_length_subaerial.name, 30_000.0, "m"),
        ]
    )
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.eruption_interval_yr.name,
        values               = [200_000.0, 300_000.0],
        units                = "yr",
        fixed                = [
            (PARAMS.extrusion_volume_factor_max.name, 0.7, "None"),
        ]
    )
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.characteristic_flow_length_subaerial.name,
        values               = [20_000.0, 30_000.0, 40_000.0],
        units                = "m",
        fixed                = [
            (PARAMS.delta_temperature_transient.name, 100.0, "deltaK"),
            (PARAMS.eruption_interval_yr.name, 200_000.0, "yr"),
            (PARAMS.extrusion_volume_factor_max.name, 0.5, "None"),
        ]
    )
    print_case_info(
        case_inputs=case_inputs, 
        case_id_max=case_id, 
        target_names=[
            PARAMS.delta_temperature_transient.name,
            PARAMS.characteristic_flow_length_subaerial.name,
            PARAMS.eruption_interval_yr.name,
            PARAMS.extrusion_volume_factor_max.name,
            PARAMS.maximum_damage_probability.name,
            PARAMS.melt_damage_distance.name,
            PARAMS.melt_damage_factor.name,
            PARAMS.dike_fluid_marker_fraction.name,
            PARAMS.density_dike_fluid.name,
            PARAMS.iuse_melt_compaction.name,
        ]
    )
    case_parameters = case_inputs[model_case_name]
    convert_case_parameters_to_standard_units!(case_parameters)
    return case_parameters
end

end # module
