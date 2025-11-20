module CaseInputTools

include("core/CaseTypes.jl")
include("core/GetInputs.jl")
include("core/CaseBuilder.jl")
include("core/CasePrinters.jl")
include("core/CaseChecker.jl")

import EarthBox.ParameterRegistry: get_parameter_dict
import .CaseChecker: convert_parameters_to_standard_units!
import .CaseTypes: CaseType
import .CaseTypes: CaseCollectionType
import .CaseTypes: CaseParameter

"""
    initialize_cases(
        base_case::CaseType,
        number_of_cases::Int = 500
    )::CaseCollectionType

Initialize a collection of model cases with a base case and a given number of cases.

# Arguments
- `base_case::CaseType`: Base case inputs.
- `number_of_cases::Int`: Number of cases to initialize.

# Returns
- `case_inputs::CaseCollectionType`: Collection of model cases.
"""
function initialize_cases(
    base_case::CaseType,
    number_of_cases::Int = 500
)::CaseCollectionType
    case_inputs = CaseCollectionType()
    for i in 1:number_of_cases
        case_name = "case$(i-1)"  # Convert to 0-based naming for consistency
        case_inputs[case_name] = deepcopy(base_case)
    end
    return case_inputs
end

""" 
    get_case_inputs_and_convert_to_standard_units(
        model_case_name::String,
        define_case_inputs::Function
    )::CaseType

Get input for active model case and convert to standard units.

# Arguments
- `model_case_name::String`
    - Name of target model case.
- `define_case_inputs::Function`
    - Function to define case inputs of type `CaseCollectionType`. 
  
# Returns
- `case_inputs_active::CaseType`
    - Dictionary of key inputs for active case with values converted to standard units.
"""
function get_case_inputs_and_convert_to_standard_units(
    model_case_name::String,
    define_case_inputs::Function
)::CaseType
    case_inputs = define_case_inputs()
    case_parameters = case_inputs[model_case_name]
    master_dict = get_parameter_dict()
    convert_parameters_to_standard_units!(case_parameters, master_dict)
    return case_parameters
end

""" 
    convert_case_parameters_to_standard_units!(
        case_parameters::CaseType
    )::Nothing

Convert case parameters to standard units.

# Arguments
- `case_parameters::CaseType`
    - Dictionary of case parameters where the key is a parameter name and the 
       value is 
       a CaseParameter struct containing the parameter value and units.

# Returns
- `nothing`
"""
function convert_case_parameters_to_standard_units!(
    case_parameters::CaseType
)::Nothing
    master_dict = get_parameter_dict()
    convert_parameters_to_standard_units!(case_parameters, master_dict)
    return nothing
end

"""
    get_value_dict(case_inputs::CaseType) -> Dict{String, Union{Float64, Int64, Bool}}

Get value dictionary for case inputs.

This function extracts the value from the parameter information collection
and builds a dictionary where keys are parameter names and values are
parameter values.

# Arguments
- `case_inputs::CaseType`: Case inputs

# Returns
- `value_dict::Dict{String, Union{Float64, Int64, Bool}}`: Dictionary where keys 
  are parameter names and values are parameter values
"""
function get_value_dict(case_inputs::CaseType)::Dict{String, Union{Float64, Int64, Bool}}
    value_dict = Dict{String, Union{Float64, Int64, Bool}}()
    for (parameter_name, collection) in case_inputs
        value_dict[parameter_name] = collection.value
    end
    return value_dict
end

end # module