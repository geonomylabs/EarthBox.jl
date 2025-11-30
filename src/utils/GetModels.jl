module GetModels

import ..PrintFuncs: print_info, print_warning

export get_models

"""
    get_models(; target_dir::Union{String, Nothing} = nothing)

Copy the EarthBox models directory to a specified target directory or the current 
working directory.

This function copies all example models from the EarthBox package to make them 
easily accessible for users to run and modify. The copy is robust and will work 
regardless of where EarthBox is installed.

# Keyword Arguments
- `target_dir::Union{String, Nothing}`: Target directory where models will be 
    copied. If `nothing` (default), copies to current working directory (pwd()).

# Returns
- `String`: Path to the copied models directory

# Examples
```julia
# Copy to current directory
EarthBox.get_models()

# Copy to a specific directory
EarthBox.get_models(target_dir = "/path/to/my/projects")
```

# Notes
- If the target models directory already exists, the function will ask for 
  confirmation before overwriting.
- The function uses `pkgdir` to robustly locate the EarthBox package directory.
"""
function get_models(; target_dir::Union{String, Nothing} = nothing)
    # Get the EarthBox package directory robustly
    # We need to get the parent module (EarthBox) directory since this is a 
    # submodule
    earthbox_pkg_dir = pkgdir(parentmodule(@__MODULE__))
    
    if earthbox_pkg_dir === nothing
        error("Could not determine EarthBox package directory. Make sure " *
              "EarthBox is properly installed.")
    end
    
    # Construct path to models directory within the package
    models_source_dir = joinpath(earthbox_pkg_dir, "models")
    
    # Check if source models directory exists
    if !isdir(models_source_dir)
        error("Models directory not found at: $(models_source_dir). " *
              "Please check your EarthBox installation.")
    end
    
    # Determine target directory
    if target_dir === nothing
        target_dir = pwd()
    end
    
    # Create full path to destination
    models_dest_dir = joinpath(target_dir, "models")
    
    # Check if destination already exists
    if isdir(models_dest_dir)
        print_warning("Models directory already exists at: $(models_dest_dir)")
        print("Overwrite existing directory? (y/n): ")
        response = readline()
        if lowercase(strip(response)) != "y"
            print_info("Cancelled. Models directory not copied.")
            return models_dest_dir
        end
        # Remove existing directory
        rm(models_dest_dir; recursive=true, force=true)
        print_info("Removed existing models directory.")
    end
    
    # Copy the models directory
    try
        cp(models_source_dir, models_dest_dir; force=true)
        print_info("Successfully copied models directory to: $(models_dest_dir)")
        print_info("You can now explore and run the example models!")
    catch e
        error("Failed to copy models directory: $(e)")
    end
    
    return models_dest_dir
end

end # module