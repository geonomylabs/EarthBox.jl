module CasePrinters

import Printf: @sprintf
import EarthBox.PrintFuncs: PRINT_SETTINGS
import EarthBox.ParameterRegistry: get_eb_parameters
import ..CaseTypes: CaseCollectionType, CaseParameter

"""
    print_case_info(
        case_inputs::CaseCollectionType,
        case_id_max::Int,
        target_names::Vector{String}
    )::Nothing

Print case information to the console in a markdown table format.

# Arguments
- `case_inputs::CaseCollectionType`: Collection of model cases.
- `case_id_max::Int`: Maximum case ID to print.
- `target_names::Vector{String}`: List of parameter names to print.

# Returns
- `nothing`: Returns nothing.
"""
function print_case_info(;
    case_inputs::CaseCollectionType,
    case_id_max::Int,
    target_names::Vector{String}
)::Nothing
    if PRINT_SETTINGS.print_case_info == false
        return nothing
    end
    
    # Check if we have cases to print
    if isempty(case_inputs)
        println("No cases to print.")
        return nothing
    end
    
    # Check if target_names is empty
    if isempty(target_names)
        println("No parameters specified in target_names.")
        return nothing
    end
    
    # Set fixed column width
    fixed_col_width = 10
    case_col_width = 6
    
    # Collect all formatted values
    formatted_data = Vector{Vector{String}}()
    case_numbers = Int[]
    
    # Collect all data for cases
    for (case_name, case_parameters) in case_inputs
        case_num = parse(Int, case_name[5:end])
        if case_num > case_id_max
            break
        end
        
        push!(case_numbers, case_num)
        row_data = String[]
        
        for param_name in target_names
            if haskey(case_parameters, param_name)
                param = case_parameters[param_name]
                formatted_val = format_value(param.value)
                push!(row_data, formatted_val)
            else
                push!(row_data, "N/A")
            end
        end
        
        push!(formatted_data, row_data)
    end
    
    # Wrap parameter names into lines of max fixed_col_width characters
    wrapped_names = []
    max_lines = 1
    for param_name in target_names
        lines = []
        remaining = param_name
        while length(remaining) > fixed_col_width
            # Try to break at underscore or other separator if possible
            break_idx = findlast('_', remaining[1:fixed_col_width])
            if break_idx === nothing || break_idx < 3
                break_idx = fixed_col_width
            end
            push!(lines, remaining[1:break_idx])
            remaining = remaining[break_idx+1:end]
        end
        if length(remaining) > 0
            push!(lines, remaining)
        end
        push!(wrapped_names, lines)
        max_lines = max(max_lines, length(lines))
    end
    
    # Build and print multi-line header
    for line_idx in 1:max_lines
        header_parts = ["| " * rpad(line_idx == 1 ? "Case" : "", case_col_width) * " "]
        for wrapped in wrapped_names
            if line_idx <= length(wrapped)
                push!(header_parts, "| " * rpad(wrapped[line_idx], fixed_col_width) * " ")
            else
                push!(header_parts, "| " * " "^fixed_col_width * " ")
            end
        end
        println(join(header_parts, "") * "|")
    end
    
    # Build and print separator
    separator_parts = ["|" * "-"^(case_col_width + 2)]
    for i in 1:length(target_names)
        push!(separator_parts, "|" * "-"^(fixed_col_width + 2))
    end
    println(join(separator_parts, "") * "|")
    
    # Print data rows
    for (idx, case_num) in enumerate(case_numbers)
        row_parts = ["| " * rpad(string(case_num), case_col_width) * " "]
        for value in formatted_data[idx]
            # Truncate values that are too long for the fixed width
            display_value = length(value) > fixed_col_width ? value[1:fixed_col_width] : value
            push!(row_parts, "| " * rpad(display_value, fixed_col_width) * " ")
        end
        println(join(row_parts, "") * "|")
    end
    
    return nothing
end

"""
Format a parameter value for display in the markdown table.

# Arguments
- `value::Union{Float64, Int64, Bool, String, Symbol}`: Value to format

# Returns
- `String`: Formatted value string
"""
function format_value(value::Union{Float64, Int64, Bool, String, Symbol})::String
    if value isa Float64
        # Use scientific notation for very small/large numbers
        if abs(value) < 0.001 && value != 0.0
            return @sprintf("%.2e", value)
        elseif abs(value) > 10000.0
            return @sprintf("%.2e", value)
        elseif value == floor(value)
            # If it's effectively an integer, display as such
            return @sprintf("%.0f", value)
        else
            # Otherwise use reasonable precision
            return @sprintf("%.4g", value)
        end
    elseif value isa Int64
        return string(value)
    elseif value isa Bool
        return string(value)
    else  # String or Symbol
        return string(value)
    end
end

end # module 