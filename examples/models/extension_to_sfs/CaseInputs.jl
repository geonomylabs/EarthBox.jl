"""
    CaseInputs

Custom module for defining case inputs.

| Case   | temperatur | velocity_  |
|        | e_base_    | second_    |
|        | lith       | step_      |
|        |            | factor     |
|--------|------------|------------|
| 0      | 1330       | 4          |
| 1      | 1330       | 1          |
| 2      | 1330       | 1.5        |
| 3      | 1330       | 2          |
| 4      | 1330       | 3          |
| 5      | 1350       | 1          |
| 6      | 1350       | 1.5        |
| 7      | 1350       | 2          |
| 8      | 1350       | 3          |
| 9      | 1350       | 4          |

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
    # maximum full extension velocity = 2.0 cm/yr
    base_case = CaseType(
        PARAMS.temperature_base_lith.name => CaseParameter(1330.0, "C"),
        PARAMS.velocity_second_step_factor.name => CaseParameter(2.0/0.5, "None"), # 0.5 cm/yr to 2 cm/yr
    )
    # Initialize using the base case
    case_inputs = initialize_cases(base_case)
    # Variable maximum full extension velocity, base temperature = 1330 C
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = 1,
        parameter_name = PARAMS.velocity_second_step_factor.name,
        values         = [0.5/0.5, 0.75/0.5, 1.0/0.5, 1.5/0.5],
        units          = "None",
    )
    # Variable maximum full extension velocity, base temperature = 1350 C
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.velocity_second_step_factor.name,
        values         = [0.5/0.5, 0.75/0.5, 1.0/0.5, 1.5/0.5, 2.0/0.5],
        units          = "None",
        fixed          = [
            (PARAMS.temperature_base_lith.name, 1350.0, "C"),
        ]
    )
    print_case_info(
        case_inputs=case_inputs, 
        case_id_max=case_id, 
        target_names=[
            PARAMS.temperature_base_lith.name,
            PARAMS.velocity_second_step_factor.name,
        ]
    )
    case_parameters = case_inputs[model_case_name]
    convert_case_parameters_to_standard_units!(case_parameters)
    return case_parameters
end

end # module
