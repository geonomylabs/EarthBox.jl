module PathValidation

using Printf

const DEFAULT_ALLOWED_PREFIXES = ["/home", "/mnt"]

"""
    validate_path(path::String;
                  allowed_prefixes::Vector{String} = ["/home", "/mnt"],
                  resolve_symlinks::Bool = false) -> String

Validate a path to prevent directory traversal and ensure absolute paths fall
under one of the allowed prefixes. Returns the normalized path if valid, throws
an error otherwise.

# Arguments
- `path::String`: The path to validate.
- `allowed_prefixes::Vector{String}`: Allowed root prefixes for absolute paths.
   Defaults to `["/home", "/mnt"]`. Relative paths are not subject to this check.
- `resolve_symlinks::Bool`: If `true` and the path exists, resolve symlinks via
   `realpath` and re-check the prefix on the resolved path. Defaults to `false`.

# Returns
- `String`: The normalized path if valid.

# Throws
- `ArgumentError` if the path attempts directory traversal or escapes the allowed
   prefixes.
"""
function validate_path(path::String;
                       allowed_prefixes::Vector{String} = DEFAULT_ALLOWED_PREFIXES,
                       resolve_symlinks::Bool = false)::String
    normalized_path = normpath(path)

    if any(==(".."), splitpath(normalized_path))
        throw(ArgumentError("Path contains directory traversal attempt: $path"))
    end

    if isabspath(normalized_path)
        if !any(startswith(normalized_path, prefix) for prefix in allowed_prefixes)
            throw(ArgumentError("Absolute path is not in allowed prefixes: $normalized_path"))
        end

        if resolve_symlinks && ispath(normalized_path)
            resolved = realpath(normalized_path)
            if !any(startswith(resolved, prefix) for prefix in allowed_prefixes)
                throw(ArgumentError(
                    "Resolved path escapes allowed prefixes: $normalized_path -> $resolved"
                ))
            end
            return resolved
        end
    end

    return normalized_path
end

"""
    validate_directory_path(path::String;
                            allowed_prefixes::Vector{String} = ["/home", "/mnt"],
                            resolve_symlinks::Bool = false) -> String

Validate a directory path and ensure the directory exists. Forwards
`allowed_prefixes` and `resolve_symlinks` to [`validate_path`](@ref).
"""
function validate_directory_path(path::String;
                                 allowed_prefixes::Vector{String} = DEFAULT_ALLOWED_PREFIXES,
                                 resolve_symlinks::Bool = false)::String
    normalized_path = validate_path(path;
                                    allowed_prefixes = allowed_prefixes,
                                    resolve_symlinks = resolve_symlinks)

    if !isdir(normalized_path)
        throw(ArgumentError("Directory does not exist: $normalized_path"))
    end

    return normalized_path
end

"""
    validate_file_path(path::String;
                       allowed_prefixes::Vector{String} = ["/home", "/mnt"],
                       resolve_symlinks::Bool = false) -> String

Validate a file path and ensure the file exists. Forwards `allowed_prefixes`
and `resolve_symlinks` to [`validate_path`](@ref).
"""
function validate_file_path(path::String;
                            allowed_prefixes::Vector{String} = DEFAULT_ALLOWED_PREFIXES,
                            resolve_symlinks::Bool = false)::String
    normalized_path = validate_path(path;
                                    allowed_prefixes = allowed_prefixes,
                                    resolve_symlinks = resolve_symlinks)

    if !isfile(normalized_path)
        throw(ArgumentError("File does not exist: $normalized_path"))
    end

    return normalized_path
end

# ----- Output-path safety guard -----

# System directories that EarthBox should never write outputs to. The check is
# a deny-list rather than an allow-list so legitimate non-/home roots used on
# HPC/workstation setups (/scratch, /data, /opt, /tmp, /mnt, ...) keep working
# without configuration. Adding an entry here means absolute paths starting
# with that prefix are rejected; relative paths are not subject to this check.
const SYSTEM_DENY_PREFIXES = [
    "/etc",   "/usr",   "/bin",  "/sbin",
    "/boot",  "/proc",  "/sys",  "/dev",
    "/lib",   "/lib64", "/run",  "/root",
]

"""
    validate_safe_output_path(path) -> String

Validate a path that EarthBox is about to create or write to. Normalises the
path, rejects directory-traversal attempts, and refuses absolute paths that
resolve under common system directories (`/etc`, `/usr`, `/bin`, `/proc`, ...).
Relative paths and other absolute paths (`/home`, `/mnt`, `/scratch`, `/tmp`,
`/data`, `/opt`, ...) are accepted unchanged.

This is a defence against typos and obviously-misconfigured output roots, not
an exhaustive sandbox. Callers that need a strict allow-list should use
[`validate_path`](@ref) with explicit `allowed_prefixes` instead.

# Throws
- `ArgumentError` if the path is empty, contains a `..` traversal segment after
   normalisation, or resolves under a system-only prefix.
"""
function validate_safe_output_path(path::AbstractString)::String
    if isempty(path)
        throw(ArgumentError("Output path is empty."))
    end
    normalized = normpath(String(path))

    if any(==(".."), splitpath(normalized))
        throw(ArgumentError("Output path contains directory traversal: $(path)"))
    end

    if isabspath(normalized)
        for prefix in SYSTEM_DENY_PREFIXES
            if normalized == prefix || startswith(normalized, prefix * "/")
                throw(ArgumentError(
                    "Output path resolves under system directory $(prefix) and " *
                    "is not allowed for EarthBox output: $(path)"
                ))
            end
        end
    end

    return normalized
end

end # module PathValidation
