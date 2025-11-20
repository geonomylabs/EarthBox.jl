module PathValidation

using Printf

"""
    validate_path(path::String) -> String

Validate a path to prevent directory traversal attacks and ensure it's within allowed boundaries.
Returns the normalized path if valid, throws an error otherwise.

# Arguments
- `path::String`: The path to validate

# Returns
- `String`: The normalized path if valid

# Throws
- `ArgumentError` if the path contains directory traversal attempts or is invalid
"""
function validate_path(path::String)::String
    # Check for directory traversal attempts
    if contains(path, "..") || contains(path, "~")
        throw(ArgumentError("Path contains directory traversal attempt: $path"))
    end

    # Normalize the path
    normalized_path = normpath(path)

    # Check if path is absolute and starts with allowed prefixes
    if isabspath(normalized_path)
        allowed_prefixes = ["/home", "/mnt"]
        if !any(startswith(normalized_path, prefix) for prefix in allowed_prefixes)
            throw(ArgumentError("Absolute path is not in allowed prefixes: $normalized_path"))
        end
    end

    return normalized_path
end

"""
    validate_directory_path(path::String) -> String

Validate a directory path and ensure it exists.
Returns the normalized path if valid, throws an error otherwise.

# Arguments
- `path::String`: The directory path to validate

# Returns
- `String`: The normalized path if valid

# Throws
- `ArgumentError` if the path is invalid or directory doesn't exist
"""
function validate_directory_path(path::String)::String
    normalized_path = validate_path(path)
    
    if !isdir(normalized_path)
        throw(ArgumentError("Directory does not exist: $normalized_path"))
    end

    return normalized_path
end

"""
    validate_file_path(path::String) -> String

Validate a file path and ensure it exists.
Returns the normalized path if valid, throws an error otherwise.

# Arguments
- `path::String`: The file path to validate

# Returns
- `String`: The normalized path if valid

# Throws
- `ArgumentError` if the path is invalid or file doesn't exist
"""
function validate_file_path(path::String)::String
    normalized_path = validate_path(path)
    
    if !isfile(normalized_path)
        throw(ArgumentError("File does not exist: $normalized_path"))
    end

    return normalized_path
end

end # module PathValidation 