module MaterialContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")
include("dicts/DictCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import EarthBox.MaterialColorsContainer: MaterialColors
import EarthBox.PrintFuncs: print_info
import .ParameterCollection: Parameters
import .ArrayCollection
import .ArrayCollection: Arrays
import .DictCollection: Dicts

"""
    Materials <: CollectionContainer

Data structure containing parameter, array and dictionary objects associated with material properties.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for material properties
- `arrays::`[`Arrays`](@ref): Array groups containing material property values
- `dicts::`[`Dicts`](@ref): Dictionary groups for mapping material names to IDs
- `colors::`[`MaterialColors`](@ref): Struct with material color information

# Constructor
    Materials()::Materials

"""
mutable struct Materials <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
    dicts::Dicts
    colors::MaterialColors
end

function Materials()
    parameters = Parameters()
    nmats = parameters.material_description.nmats.value
    print_info("Maximum number of allowed materials (nmats): $nmats")
    arrays = Arrays(nmats)
    dicts = Dicts(
        Dict{String, Int64}(),
        Dict{String, Vector{Int64}}()
    )
    colors = MaterialColors()
    return Materials(parameters, arrays, dicts, colors)
end

function override_plasticity_in_zones_with_strain_weakening!(
    state::Materials,
    strain_final::Union{Float64, Nothing},
    friction_angle_initial::Union{Float64, Nothing},
    friction_angle_final::Union{Float64, Nothing},
    material_ids_zones_with_strain_weakening::Vector{Int16}
)
    for material_id in material_ids_zones_with_strain_weakening
        ArrayCollection.set_final_strain(
            state.arrays, material_id, strain_final)
        ArrayCollection.set_friction_angle_initial(
            state.arrays, material_id, friction_angle_initial)
        ArrayCollection.set_friction_angle_final(
            state.arrays, material_id, friction_angle_final)
    end
end

function override_plasticity_in_strong_zones!(
    state::Materials,
    strain_final::Union{Float64, Nothing},
    friction_angle_initial::Union{Float64, Nothing},
    friction_angle_final::Union{Float64, Nothing},
    material_ids_strong_zones::Vector{Int16}
)
    for material_id in material_ids_strong_zones
        ArrayCollection.set_final_strain(
            state.arrays, material_id, strain_final)
        ArrayCollection.set_friction_angle_initial(
            state.arrays, material_id, friction_angle_initial)
        ArrayCollection.set_friction_angle_final(
            state.arrays, material_id, friction_angle_final)
    end
end

function override_plasticity!(
    state::Materials,
    strain_final::Union{Float64, Nothing},
    friction_angle_initial::Union{Float64, Nothing},
    friction_angle_final::Union{Float64, Nothing},
    material_ids_target::Vector{Int16}
)
    for material_id in material_ids_target
        ArrayCollection.set_final_strain(
            state.arrays, material_id, strain_final)
        ArrayCollection.set_friction_angle_initial(
            state.arrays, material_id, friction_angle_initial)
        ArrayCollection.set_friction_angle_final(
            state.arrays, material_id, friction_angle_final)
    end
end

function override_pre_exponential_factors_in_asthenosphere!(
    state::Materials,
    scale_factor_for_dislocation_creep::Union{Float64, Nothing},
    scale_factor_for_diffusion_creep::Union{Float64, Nothing},
    material_ids_asthenosphere::Vector{Int16}
)
    for material_id in material_ids_asthenosphere
        ArrayCollection.scale_pre_exponential_factor_dislocation_creep(
            state.arrays, material_id, scale_factor_for_dislocation_creep)
        ArrayCollection.scale_pre_exponential_factor_diffusion_creep(
            state.arrays, material_id, scale_factor_for_diffusion_creep)
    end
end

function override_pre_exponential_factors_in_continental_crust!(
    state::Materials,
    scale_factor_for_dislocation_creep::Union{Float64, Nothing},
    material_ids_continental_crust::Vector{Int16}
)
    for material_id in material_ids_continental_crust
        ArrayCollection.scale_pre_exponential_factor_dislocation_creep(
            state.arrays, material_id, scale_factor_for_dislocation_creep)
    end
end

function override_pre_exponential_factors_dislocation_creep!(
    state::Materials,
    scale_factor_for_dislocation_creep::Union{Float64, Nothing},
    material_ids_target::Vector{Int16}
)
    for material_id in material_ids_target
        ArrayCollection.scale_pre_exponential_factor_dislocation_creep(
            state.arrays, material_id, scale_factor_for_dislocation_creep)
    end
end

function override_radiogenic_heat_production_in_continental_crust!(
    state::Materials,
    radiogenic_heat_production_felsic_crust::Union{Float64, Nothing},
    radiogenic_heat_production_mafic_crust::Union{Float64, Nothing},
    material_ids_felsic_continental_crust::Vector{Int16},
    material_ids_mafic_continental_crust::Vector{Int16}
)
    for material_id in material_ids_felsic_continental_crust
        ArrayCollection.set_radiogenic_heat_production(
            state.arrays, material_id, radiogenic_heat_production_felsic_crust)
    end

    for material_id in material_ids_mafic_continental_crust
        ArrayCollection.set_radiogenic_heat_production(
            state.arrays, material_id, radiogenic_heat_production_mafic_crust)
    end
end

function override_latent_heat!(
    state::Materials,
    latent_heat_new::Union{Float64, Nothing},
    material_ids_target::Vector{Int16}
)
    for material_id in material_ids_target
        ArrayCollection.set_latent_heat(
            state.arrays, material_id, latent_heat_new)
    end
end

function override_solidus!(
    state::Materials,
    solidus_model_new::Union{String, Nothing},
    material_ids_target::Vector{Int16}
)
    for material_id in material_ids_target
        ArrayCollection.set_solidus_model(
            state.arrays, material_id, solidus_model_new)
    end
end

function override_liquidus!(
    state::Materials,
    liquidus_model_new::Union{String, Nothing},
    material_ids_target::Vector{Int16}
)
    for material_id in material_ids_target
        ArrayCollection.set_liquidus_model(
            state.arrays, material_id, liquidus_model_new)
    end
end

end # module 