module MaterialCollectionManager

import EarthBox.PrintFuncs: print_info, print_warning

export get_material_collection

"""
    MaterialCollection

# Fields
- `materials::NamedTuple`:
    - Named tuple of material names associated with the material collection.
- `path::String`:
    - Path to the yaml formatted material collection file with groups of 
      materials and their properties.
"""
mutable struct MaterialCollection
    materials::NamedTuple
    path::String
end

"""
    get_material_collection(
        collection::MaterialCollection; 
        target_dir::Union{String, Nothing} = nothing
    )

Copy a material library collection .yml file from a MaterialCollection to a 
specified target directory or the current working directory.

This function copies material library files to make them easily accessible for 
users to view, modify, and use in their models. The copy is robust and will work 
regardless of where EarthBox is installed.

# Arguments
- `collection::MaterialCollection`: A MaterialCollection object whose `path` 
    attribute points to the .yml file to be copied.

# Keyword Arguments
- `target_dir::Union{String, Nothing}`: Target directory where the file will be 
    copied. If `nothing` (default), copies to current working directory (pwd()).

# Returns
- `String`: Full path to the copied file

# Examples
```julia
using EarthBox
# Get a material collection and copy its file to current directory
mat_lib = MaterialLibrary()
get_material_collection(mat_lib.lithospheric_deformation.lithospheric_deformation_eb1)

# Copy to specific directory
get_material_collection(
    mat_lib.sandbox,
    target_dir = "/path/to/my/project"
)
```

# Notes
- If the target file already exists, the function will ask for confirmation 
  before overwriting.
- The original filename is preserved.
"""
function get_material_collection(
    collection::MaterialCollection; 
    target_dir::Union{String, Nothing} = nothing
)
    # Get the source path from the MaterialCollection
    source_path = collection.path
    
    # Validate that the source path exists and is a file
    if !isfile(source_path)
        error("Source file not found at: $(source_path). " *
              "Please check the MaterialCollection path.")
    end
    
    # Validate that the file has .yml extension
    if !endswith(source_path, ".yml")
        error("Source file must have .yml extension. Got: $(source_path)")
    end
    
    # Determine target directory
    if target_dir === nothing
        target_dir = pwd()
    end
    
    # Create target directory if it doesn't exist
    if !isdir(target_dir)
        print_info("Creating target directory: $(target_dir)")
        mkpath(target_dir)
    end
    
    # Use original filename
    target_filename = basename(source_path)
    
    # Create full path to destination
    target_path = joinpath(target_dir, target_filename)
    
    # Check if destination already exists
    if isfile(target_path)
        print_warning("File already exists at: $(target_path)")
        print("Overwrite existing file? (y/n): ")
        response = readline()
        if lowercase(strip(response)) != "y"
            print_info("Cancelled. File not copied.")
            return target_path
        end
    end
    
    # Copy the file
    try
        cp(source_path, target_path; force=true)
        print_info("Successfully copied material library file to: $(target_path)")
    catch e
        error("Failed to copy file: $(e)")
    end
    
    return target_path
end

end