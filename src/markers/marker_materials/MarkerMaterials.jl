module MarkerMaterials

include("utils/Registry.jl")
include("utils/GetProperty.jl")
include("utils/MaterialGroupIDs.jl")
include("options/Options.jl")
include("init_manager/InitManager.jl")
include("boundary_friction/MarkerBoundaryFriction.jl")
include("stress_limits/MarkerStressLimits.jl")
include("viscous_strain_softening/MarkerViscousStrainSoftening.jl")
include("materials/MaterialsContainer.jl")
include("get_ids/GetMaterialIDs.jl")
include("override/Override.jl")
include("utils/MaterialOverride.jl")

import EarthBox.PrintFuncs: print_info
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: InputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox: InitializationTools
import EarthBox.ParameterGroupTools: set_group_parameters!
import EarthBox.UseOptionTools: set_use_option!
import EarthBox: OptionTools
import .MaterialsContainer
import .Options: get_options
import .Options: option_ids
import .Options: option_names
import .InitManager
import .MarkerBoundaryFriction
import .MarkerStressLimits
import .MarkerViscousStrainSoftening

export initialize!

const MAT_OPTIONS = get_options()

const PDATA = get_eb_parameters()

Base.@kwdef struct PathKeyNames
    materials_library_file::String = "materials_library_file"
    materials_input_file::String = "materials_input_file"
end

struct ValidInputNames
    viscosity_min::Symbol
    viscosity_max::Symbol
    iuse_fluid_pressure_for_yield::Symbol
    plastic_healing_rate::Symbol
end

function make_material_models_string()::String
    material_model_string = ""
    
    for (option_id, option_state) in MAT_OPTIONS
        option_name = Symbol(option_state.option_name)
        required_geometries = option_state.required_geometries
        list_rgeometries = join(["    - $geometry" for geometry in required_geometries], "\n")
        
        material_model_string *= """
## $(option_state.option_name)
- `material_model` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
- **Description**: $(option_state.description)
- **Required Domains**: $(join(option_state.required_domains, ", "))
- **Required Geometries**:
$list_rgeometries
"""
    end
    
    return material_model_string
end

"""
    initialize!(
        model::ModelData;
        material_model::Union{Int, String, Symbol, Nothing}=nothing,
        paths::Union{Dict{String, String}, Nothing}=nothing,
        materials_input_dict::Union{MaterialsDictType, Nothing}=nothing,
        kwargs...
    )::Nothing

Initialize marker materials.

# Arguments
- `model::ModelData`: The model data object containing model parameters and arrays.
- `material_model::Union{Int, String, Symbol, Nothing}`: Controls the type of 
    material model used to define material properties for each marker. See the **Material Models** 
    section below for information on available material models. The material model is stored in the 
    model data container as an integer ID (`itype_mat`) and a corresponding string name 
    (`stype_mat`). If `material_model` is nothing, the current material model defined in the model 
    data container will be used. Material model parameters can be accessed from the model data 
    container as follows:
    - `itype_mat = model.materials.parameters.material_description.itype_mat.value`
    - `stype_mat = model.materials.parameters.material_description.stype_mat.value`
- `paths::Union{Dict{String, String}, Nothing}`: The paths to the materials 
    library and material input files. Dictionary keys and value descriptions are as follows:
    - `"materials_library_file"`: Path to a yaml formatted material collection library file with
         material properties for each material. The material_library_file path is required. See 
         [`MaterialLibrary`](@ref EarthBox.MaterialLibraryCollection.MaterialLibrary) for more 
         information about using material collection packaged with EarthBox. See 
         [Material Collection Files](@ref) for more information about creating a custom material 
         collection file.
    - `"materials_input_file"`: Path to a material input file with yaml 
         formatted material input definitions for each material used in the model including 
         material ID, material name from the material library, material type, material 
         domain, and red, green and blue color values. This key is optional. If not provided,
         the `materials_input_dict` must be provided. See [Material Input Files](@ref) for more information 
         about material input files.
- `materials_input_dict::Union{MaterialsDictType, Nothing}`: The materials input dictionary
    containing material input definitions for each material used in the model including 
    material ID, material name from the material library, material type, material 
    domain, and red, green and blue color values. This key is optional. If not provided,
    the `materials_input_file` must be provided. See [Material Input Dictionaries](@ref) for more information 
    about material input dictionaries.

# Keyword Arguments
- `viscosity_min::Float64`: 
    - $(PDATA.viscosity_min.description)
- `viscosity_max::Float64`:
    - $(PDATA.viscosity_max.description)
- `iuse_fluid_pressure_for_yield::Int`:
    - $(PDATA.iuse_fluid_pressure_for_yield.description)
- `plastic_healing_rate::Float64`:
    - $(PDATA.plastic_healing_rate.description)

---
# Material Models
---
This section provides a detailed description of each material model available in EarthBox
including alternative `material_model` values that can be used to initialize the material model,
a description of the material model, and the required domains and geometries for each material model.
For more information about required domains see 
[MaterialDomainsRegistry](@ref EarthBox.Markers.MarkerMaterials.Registry.MaterialDomainsRegistry).
See the links for required geometries for more information about setting required geometric 
parameters for each material model.

If the domain "MantleLithosphere" is defined in the material model then the domains 
"UpperMantleLithosphere", "MiddleMantleLithosphere" and "LowerMantleLithosphere"
will have properties assigned to them automatically from the "MantleLithosphere" domain.

$(make_material_models_string())

"""
function initialize!(
    model::ModelData;
    material_model::Union{Int, String, Symbol, Nothing}=nothing,
    paths::Union{Dict{String, String}, Nothing}=nothing,
    materials_input_dict::Union{MaterialsDictType, Nothing}=nothing,
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    option_id = InitializationTools.update_option_id_using_input_option_name(
        get_options(), material_model, model,
        get_option_id_from_model, update_option_id
        )
    check_paths!(paths)
    load_materials(option_id, model, paths, materials_input_dict)
    option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    InitManager.initialize!(model, Val(option_name))
    return nothing
end

function load_materials(
    option_id::Int,
    model::ModelData,
    paths::Dict{String, String},
    materials_input_dict::Union{MaterialsDictType, Nothing}
)::Nothing
    materials = MaterialsContainer.Materials()
    MaterialsContainer.set_material_option!(materials, option_id, get_options()[option_id])
    MaterialsContainer.load!(
        materials,
        model                      = model,
        material_library_file_path = paths["materials_library_file"],
        material_model_file_path   = paths["materials_input_file"],
        materials_input_dict       = materials_input_dict
        )
    return nothing
end

"""
    check_paths!(paths::Union{Dict{String, String}, Nothing})::Nothing

Check the paths dictionary for the materials library and material input files.

If the `"materials_input_file"` key is not provided, it is set to an empty string.

If the `"materials_library_file"` key is not provided, an error is thrown since
the materials library file is required.

# Arguments
- `paths::Union{Dict{String, String}, Nothing}`: The paths dictionary.

# Returns
- `Nothing`: Nothing.
"""
function check_paths!(paths::Union{Dict{String, String}, Nothing})::Nothing
    if paths === nothing
        paths = Dict{String, String}()
    end
    if !haskey(paths, "materials_library_file")
        throw(ArgumentError(
            "The materials_library_file path is not defined in the paths dictionary. Adjust inputs."
        ))
    end
    if !haskey(paths, "materials_input_file")
        paths["materials_input_file"] = ""
    end
    return nothing
end

function get_option_id_from_model(model::ModelData)::Int
    return model.materials.parameters.material_description.itype_mat.value
end

function get_stype_from_model(model::ModelData)::String
    return model.materials.parameters.material_description.stype_mat.value
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.materials.parameters.material_description.itype_mat.value = option_id
    return nothing
end

function print_option(model::ModelData)::Nothing
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Marker Materials Option")
    return nothing
end

end # module 