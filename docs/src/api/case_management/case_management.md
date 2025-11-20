# Case Management

EarthBox provides a case definition system that allows the user to build a collection
of model cases and run simulations based on a standard case name like `case0`, `case1` etc.

The recommended workflow is to build a module called `CaseInputs` in a file located
within the model directory called `CaseInputs.jl` where a dictionary of model cases
is generated. This module can then be included in the a main model file `Model.jl`
providing access to model parameters based on case name.

## Case Management Data Types


```@docs
EarthBox.CaseParameter
```

```@docs
EarthBox.CaseType
```

```@docs
EarthBox.CaseCollectionType
```

## Case Builder Tools

```@docs
EarthBox.initialize_cases
```

```@docs
EarthBox.define_case_group!
```

```@docs
EarthBox.convert_case_parameters_to_standard_units!
```

```@docs
EarthBox.print_case_info
```

# Example

The following example is of a module that would be defined in a file called 
`CaseInputs.jl` and includes a function called `define_case_inputs` that returns
a dictionary of model parameters based on a specified case name. This particular 
example involves building a range of scenarios that explore the melt damage model.

```julia
using EarthBox
PRINT_SETTINGS.print_case_info = true;

module CaseInputs

using EarthBox

const PARAMS = get_eb_parameters()

function define_case_parameters(;model_case_name::String="case0")::CaseType
    # Define the base case
    base_case = CaseType(
        # ThermalBottomTransientBoundaryCondition
        PARAMS.delta_temperature_transient.name => CaseParameter(100.0, "deltaK"),
        # MeltDamageModel
        PARAMS.iuse_melt_damage.name     => CaseParameter(1, "None"),
        PARAMS.melt_damage_distance.name => CaseParameter(2.5, "km"),
        PARAMS.melt_damage_factor.name   => CaseParameter(10.0, "None"),
        # Extrusion Model
        PARAMS.extrusion_volume_factor.name               => CaseParameter(0.06, "None"),
        PARAMS.extrusion_volume_factor_max.name           => CaseParameter(0.5, "None"),
        PARAMS.characteristic_flow_length_subaerial.name  => CaseParameter(20.0, "km"),
        PARAMS.eruption_interval_yr.name                  => CaseParameter(50_000.0, "yr")
    )
    # Initialize all model cases using the base case (default is 500 cases)
    case_inputs = initialize_cases(base_case)
    # Update case_inputs with 4 cases involving variable delta temperature transient
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
    # Update case_inputs with 4 cases with variable melt damage distance
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.melt_damage_distance.name,
        values         = [0.625, 1.25, 5.0],
        units          = "km"
    )
    # Update case_inputs with 4 cases with variable eruption duration
    case_id = define_case_group!(
        case_inputs,
        case_id_ini    = case_id + 1,
        parameter_name = PARAMS.eruption_interval_yr.name,
        values         = [12_500.0, 25_000.0, 100_000.0],
        units          = "yr"
    )
    # Variable maximum extrusion volume factor @ eruption duration = 50_000 yrs
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
    # Update case_inputs with variable maximum extrusion volume factor @ eruption duration = 75_000 yrs
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.extrusion_volume_factor_max.name,
        values               = [0.1, 0.3, 0.5, 0.7],
        units                = "None",
        fixed_parameter_name = PARAMS.eruption_interval_yr.name,
        fixed_value          = 75_000.0,
        fixed_units          = "yr"
    )
    # Update case_inputs with variable maximum extrusion volume factor @ eruption duration = 100_000 yrs
    case_id = define_case_group!(
        case_inputs,
        case_id_ini          = case_id + 1,
        parameter_name       = PARAMS.extrusion_volume_factor_max.name,
        values               = [0.1, 0.2, 0.3, 0.7],
        units                = "None",
        fixed_parameter_name = PARAMS.eruption_interval_yr.name,
        fixed_value          = 100_000.0,
        fixed_units          = "yr"
    )
    # Update case_inputs with variable maximum extrusion volume factor @ damage factor = 1
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
    # Only the specified target names will show up in the markdown table
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
            PARAMS.eruption_interval_yr.name
        ]
    )
    # Select the active case
    case_parameters = case_inputs[model_case_name]
    # Convert to standard units
    convert_case_parameters_to_standard_units!(case_parameters)
    return case_parameters
end

end # module

case_parameters = CaseInputs.define_case_parameters(model_case_name="case3")

```

```
| Case   | delta_     | extrusion_ | extrusion_ | iuse_melt_ | melt_      | melt_      | characteri | eruption_  |
|        | temperatur | volume_    | volume_    | damage     | damage_    | damage_    | stic_flow_ | interval_  |
|        | e_transien | factor     | factor_max |            | factor     | distance   | length_    | yr         |
|        | t          |            |            |            |            |            | subaerial  |            |
|--------|------------|------------|------------|------------|------------|------------|------------|------------|
| 0      | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 1      | 50         | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 2      | 75         | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 3      | 125        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 4      | 150        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 5      | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 6      | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 7      | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 8      | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 9      | 100        | 0.06       | 0.5        | 1          | 1          | 2.5        | 20         | 5.00e+04   |
| 10     | 100        | 0.06       | 0.5        | 1          | 1.25       | 2.5        | 20         | 5.00e+04   |
| 11     | 100        | 0.06       | 0.5        | 1          | 2.5        | 2.5        | 20         | 5.00e+04   |
| 12     | 100        | 0.06       | 0.5        | 1          | 5          | 2.5        | 20         | 5.00e+04   |
| 13     | 100        | 0.06       | 0.5        | 1          | 10         | 0.625      | 20         | 5.00e+04   |
| 14     | 100        | 0.06       | 0.5        | 1          | 10         | 1.25       | 20         | 5.00e+04   |
| 15     | 100        | 0.06       | 0.5        | 1          | 10         | 5          | 20         | 5.00e+04   |
| 16     | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 1.25e+04   |
| 17     | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 2.50e+04   |
| 18     | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 1.00e+05   |
| 19     | 100        | 0.06       | 0.1        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 20     | 100        | 0.06       | 0.2        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 21     | 100        | 0.06       | 0.3        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 22     | 100        | 0.06       | 0.4        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 23     | 100        | 0.06       | 0.6        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 24     | 100        | 0.06       | 0.7        | 1          | 10         | 2.5        | 20         | 5.00e+04   |
| 25     | 100        | 0.06       | 0.1        | 1          | 10         | 2.5        | 20         | 7.50e+04   |
| 26     | 100        | 0.06       | 0.3        | 1          | 10         | 2.5        | 20         | 7.50e+04   |
| 27     | 100        | 0.06       | 0.5        | 1          | 10         | 2.5        | 20         | 7.50e+04   |
| 28     | 100        | 0.06       | 0.7        | 1          | 10         | 2.5        | 20         | 7.50e+04   |
| 29     | 100        | 0.06       | 0.1        | 1          | 10         | 2.5        | 20         | 1.00e+05   |
| 30     | 100        | 0.06       | 0.2        | 1          | 10         | 2.5        | 20         | 1.00e+05   |
| 31     | 100        | 0.06       | 0.3        | 1          | 10         | 2.5        | 20         | 1.00e+05   |
| 32     | 100        | 0.06       | 0.7        | 1          | 10         | 2.5        | 20         | 1.00e+05   |
| 33     | 100        | 0.06       | 0.1        | 1          | 1          | 2.5        | 20         | 5.00e+04   |
| 34     | 100        | 0.06       | 0.2        | 1          | 1          | 2.5        | 20         | 5.00e+04   |
| 35     | 100        | 0.06       | 0.3        | 1          | 1          | 2.5        | 20         | 5.00e+04   |
| 36     | 100        | 0.06       | 0.4        | 1          | 1          | 2.5        | 20         | 5.00e+04   |
```

The function `define_case_inputs` can be accessed from a local `Model.jl` file
using the following julia code:

```julia
# Model.jl
using EarthBox
include("CaseInputs.jl")
import .CaseInputs: define_case_parameters

const PARAMS = get_eb_parameters()

# Get the melt_damage_distance parameter for case3
case_parameters = define_case_parameters(model_case_name="case3")
melt_damage_distance = case_parameters[PARAMS.melt_damage_distance].value
```

# Material Override Using CaseType

Some material properties defined in the material collection can be overridden
using the API and `case_parameters` dictionary with type [EarthBox.CaseType](@ref).
Standard earthbox parameters are defined from `get_eb_parameters()` and values
for these parameters are set in the `case_parameters` dictionary. Material
overriding is based on material types and domains (see [Material Types](@ref) and
[Material Domains](@ref)). The docstring below provides a detailed descriptions
of parameters and associated material types and domains.

```@docs
Markers.MarkerMaterials.MaterialOverride.override_material_properties!
```