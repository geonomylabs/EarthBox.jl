module PrintFuncs

import LinearAlgebra
import LinearAlgebra: norm
using Printf

Base.@kwdef mutable struct Settings
    print_performance::Bool = false
    print_info::Bool = true
    print_debug::Bool = false
    print_warning::Bool = false
    print_earthbox_error::Bool = false
    print_case_info::Bool = false
    print_flow_info::Bool = false
    print_melt_extraction_info::Bool = false
end

const PRINT_SETTINGS = Settings()

const PREFIX_PERFORMANCE    = "++"
const PREFIX_DEBUG          = "** DEBUG:"
const PREFIX_WARNING        = "!!! WARNING !!!"
const PREFIX_EARTHBOX_ERROR = "!!! EarthBox Error !!!"

const INDENT = Dict(1 => "", 2 => "   ", 3 => "      ")
const BULLETS = Dict(1 => ">>", 2 => "+", 3 => "*")

"""
    @timeit_memit(description, expr)

Macro for timing and memory output.

# Arguments
- `description`: Description of the section of code being timed.
- `expr`: The expression to time and measure memory for.

# Example
```julia
@timeit_memit "Building matrix" begin
    matrix = zeros(Float64, 1000, 1000)
    # code being timed
end
```

"""
macro timeit_memit(description, expr)
    # if no level is specified assume it is called in a loop
    return _timeit_memit(description, 2, expr)
end

macro timeit_memit(description, level, expr)
    return _timeit_memit(description, level, expr)
end

function _timeit_memit(description, level, expr)
    return quote
        # Get memory usage before
        local description_str = $(esc(description))
        local mem_mb = Sys.maxrss() / 1024 / 1024  # Convert to MB
        local level_val = $(esc(level))
        # Time the function and capture its result
        local result = nothing
        local time1 = time()
        result = $(esc(expr))  # Store the result
        local time2 = time()
        local delta_time = time2 - time1
        print_performance(description_str, delta_time, mem_mb, level_val)
        result  # Explicitly return the result
    end
end


"""
# Arguments
- `description::String`: Description of the section of code being timed.
- `delta_time::Float64`: Time taken to run the section of code (seconds).
- `mem_mb::Float64`: Memory used by the section of code (MB).
"""
function print_performance(
    description::String, 
    delta_time::Float64, 
    mem_mb::Float64,
    level::Int
)
    @assert level in 1:3 "Level must be between 1 and 3"
    if PRINT_SETTINGS.print_performance
        PREFIX = INDENT[level] * PREFIX_PERFORMANCE
        @printf("%s %-80s : cpu(s) : %10.3f : MB : %6.2f\n", PREFIX, description, delta_time, mem_mb)
    end
    flush(stdout)
end

function print_melt_extraction_info(str_info::String; level::Int=1)
    if PRINT_SETTINGS.print_melt_extraction_info
        print_info(str_info, level=level)
    end
end

function print_flow_info(str_info::String; level::Int=1)
    if PRINT_SETTINGS.print_flow_info
        print_info(str_info, level=level)
    end
end

function print_info(str_info::String; level::Int=1)
    @assert level in 1:3 "Level must be between 1 and 3"
    if PRINT_SETTINGS.print_info
        println("$(INDENT[level])$(BULLETS[level]) $str_info")
    end
    flush(stdout)
end

function print_info(
    str_description::String, 
    val_list::Vector{String},
    level::Int=1
)
    @assert level in 1:3 "Level must be between 1 and 3"
    if PRINT_SETTINGS.print_info
        str_out = join(val_list, " : ")
        println("$(INDENT[level])$(BULLETS[level]) $str_description : $str_out")
    end
    flush(stdout)
end

function print_warning(str_warning::String; level::Int=1)
    @assert level in 1:3 "Level must be between 1 and 3"
    if PRINT_SETTINGS.print_warning
        println("$(INDENT[level])$PREFIX_WARNING $str_warning")
    end
    flush(stdout)
end

function print_earthbox_error(str_error::String; level::Int=1)
    @assert level in 1:3 "Level must be between 1 and 3"
    if PRINT_SETTINGS.print_earthbox_error
        println("$(INDENT[level])$PREFIX_EARTHBOX_ERROR $str_error")
    end
    flush(stdout)
end

function print_debug(str_debug::String; level::Int=1)
    @assert level in 1:3 "Level must be between 1 and 3"
    if PRINT_SETTINGS.print_debug
        println("$(INDENT[level])$PREFIX_DEBUG $str_debug")
    end
    flush(stdout)
end

function print_unit_conversion_info(
    parameter_name::String, 
    unit_start::String, 
    unit_end::String;
    level::Int=1
)
    @assert level in 1:3 "Level must be between 1 and 3"
    if PRINT_SETTINGS.print_info
        print_info("$parameter_name : $unit_start -> $unit_end", level=level)
    end
    flush(stdout)
end

function print_unit_conversion_info_and_values(
    parameter_name::String, 
    value_start::Float64,
    value_end::Float64,
    unit_start::String, 
    unit_end::String;
    level::Int=1
)
    @assert level in 1:3 "Level must be between 1 and 3"
    if PRINT_SETTINGS.print_info
        print_info("$parameter_name : $value_start $unit_start -> $value_end $unit_end", level=level)
    end
    flush(stdout)
end

function print_option_category_info(category_name::String)
    if PRINT_SETTINGS.print_info
        print_info("Option Category: $category_name", level=1)
    end
    flush(stdout)
end

function print_option_id(option_id::Int)
    print_info("Option ID: $option_id", level=2)
end

function print_option_info(option_name::String, description::String)
    if PRINT_SETTINGS.print_info
        print_info("Option Name: $option_name", level=2)
        print_info("Option Description: $description", level=2)
    end
    flush(stdout)
end

function print_list(list_name::String, data_list::Vector{Any})
    println()
    println(list_name)
    for (i, val) in enumerate(data_list)
        println("$i : $val")
    end
    println()
end

function print_list(list_name::String, data_list::Vector{String})
    println()
    println(list_name)
    for (i, val) in enumerate(data_list)
        println("$i : $val")
    end
    println()
end

function print_solution_vector_statistics(S::Vector{Float64})
    n = length(S)
    min = LinearAlgebra.minimum(S)
    max = LinearAlgebra.maximum(S)
    avg = sum(S) / n
    #print_info("Solution vector statistics:", level=1)
    #print_info("Length: $n", level=2)
    #print_info("Mean: $(round(avg, digits=6))", level=2)
end

end # module 