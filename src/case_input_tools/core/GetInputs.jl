module GetInputs

""" Get istart and iend for plotting from command line arguments.

# Arguments
- `istart::Int`: Default starting index for plotting (default: 0)
- `iend::Union{Int, Nothing}`: Default ending index for plotting (default: nothing)

# Returns
- `Tuple{Int, Union{Int, Nothing}}`: A tuple containing:
    - istart: The index from which to start plotting
    - iend: The index at which to end plotting (or nothing if not specified)

# Note
This function reads from ARGS[2] and ARGS[3] for istart and iend respectively.
"""
function get_istart_iend(;
    istart::Int=0,
    iend::Union{Int, Nothing}=nothing
)::Tuple{Int, Union{Int, Nothing}}
    if length(ARGS) > 1
        try
            istart = parse(Int, ARGS[2])
        catch
            # If parsing fails, keep the default value
        end
    end
    
    if length(ARGS) > 3
        try
            iend = parse(Int, ARGS[3])
        catch
            # If parsing fails, keep the default value
        end
    end
    
    return (istart, iend)
end

export get_istart_iend

end # module 