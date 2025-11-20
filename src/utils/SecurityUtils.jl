module SecurityUtils

import EarthBox.PrintFuncs: print_warning
using Printf
using Base.Filesystem
import Base: Union, Regex

"""
    check_file_permissions(path::String, required_permissions::String) -> Bool

Check if a file or directory has the required permissions.
Returns true if permissions are sufficient, false otherwise.

# Arguments
- `path::String`: Path to the file or directory to check
- `required_permissions::String`: String of required permissions (e.g., "rwx" for read, write, execute)

# Returns
- `Bool`: true if permissions are sufficient, false otherwise
"""
function check_file_permissions(path::String, required_permissions::String)::Bool
    for perm in required_permissions
        if perm == 'r' && !isreadable(path)
            print_warning("Missing read permission for: $path")
            return false
        elseif perm == 'w' && !iswritable(path)
            print_warning("Missing write permission for: $path")
            return false
        elseif perm == 'x' && !isexecutable(path)
            print_warning("Missing execute permission for: $path")
            return false
        end
    end
    
    return true
end

"""
    check_environment_variables(required_vars::Vector{String}) -> Dict{String, String}

Check if required environment variables are set and not empty.
Returns a dictionary of found variables and their values.

# Arguments
- `required_vars::Vector{String}`: List of required environment variable names

# Returns
- `Dict{String, String}`: Dictionary of found variables and their values

"""
function check_environment_variables(required_vars::Vector{String})::Dict{String, String}
    found_vars = Dict{String, String}()
    missing_vars = String[]
    
    for var in required_vars
        value = get(ENV, var, "")
        if isempty(value)
            push!(missing_vars, var)
        else
            found_vars[var] = value
        end
    end
    
    if !isempty(missing_vars)
        print_warning(
            "Missing environment variables (This may be OK if you are using Julia Binaries): $(join(missing_vars, ", "))", level=1
            )
    end
    
    return found_vars
end

end # module SecurityUtils 