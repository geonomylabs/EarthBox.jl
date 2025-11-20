module MaterialOverride

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: MaterialContainer
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.CaseInputTools.CaseTypes: CaseType
import ..MaterialGroupIDs: get_ids_for_plastic_materials_with_strain_weakening
import ..MaterialGroupIDs: get_ids_for_lithospheric_strong_zones
import ..MaterialGroupIDs: get_ids_for_felsic_and_mafic_continental_crust
import ..MaterialGroupIDs: get_ids_for_asthenospheric_mantle
import ..MaterialGroupIDs: get_ids_for_all_mantle_rocks
import ..MaterialGroupIDs: get_ids_for_solidified_basalt
import ..MaterialGroupIDs: get_ids_for_oceanic_crust
import ..MaterialGroupIDs: get_ids_for_sediment

const keys = get_eb_parameters()

const PARAMS = get_eb_parameters()

"""
    override_material_properties!(model::ModelData, case_inputs_active::CaseType)

Override material properties for model materials using a case type.

# Arguments
- `model::ModelData`: The model data object that contains the materials to override
- `case_inputs_active::CaseType`: The case inputs active with keys equal to
    material properties to override.

# Material Domains, Genetic Types, and Parameters for Overriding

## Plasticity in Solidified Basalt
- type_names: "SolidifiedBasalt"
- Parameters:
    - `$(PARAMS.strain_final_solidified_basalt.name)`
    - `$(PARAMS.friction_angle_initial_solidified_basalt.name)`
    - `$(PARAMS.friction_angle_final_solidified_basalt.name)`
## Plasticity in Oceanic Crust
- type_names:
    - "SolidifiedGabbro"
    - "SolidifiedGabbroPartiallyMolten"
    - "SolidifiedLayeredGabbro"
    - "SolidifiedLayeredGabbroPartiallyMolten"
- Parameters:
    - `$(PARAMS.strain_final_oceanic_crust.name)`
    - `$(PARAMS.friction_angle_initial_oceanic_crust.name)`
    - `$(PARAMS.friction_angle_final_oceanic_crust.name)`
## Plasticity in Sediment
- type_names: "Sediment"
- Parameters:
    - `$(PARAMS.strain_final_sediment.name)`
    - `$(PARAMS.friction_angle_initial_sediment.name)`
    - `$(PARAMS.friction_angle_final_sediment.name)`
## Plasticity in Zones with Strain Weakening
- type_names:
    - "FelsicContinentalCrustFertile"
    - "FelsicContinentalCrustPartiallyMolten"
    - "FelsicContinentalCrustRefactory"
    - "MaficContinentalCrustFertile"
    - "MaficContinentalCrustPartiallyMolten"
    - "MaficContinentalCrustRefactory"
    - "UltramaficMantleFertile"
    - "UltramaficMantlePartiallyMolten"
    - "UltramaficMantleRefactory"
    - "SolidifiedGabbro"
    - "SolidifiedGabbroPartiallyMolten"
    - "SolidifiedLayeredGabbro"
    - "SolidifiedLayeredGabbroPartiallyMolten"
    - "SolidifiedBasalt"
    - "Sediment"
- Parameters:
    - `$(PARAMS.strain_final.name)`
    - `$(PARAMS.friction_angle_initial.name)`
    - `$(PARAMS.friction_angle_final.name)`
### Plasticity in Lithospheric Strong Zones
- domain_names:
    - "LithosphericMantleStrongZone"
    - "UpperContinentalCrustStrongZone"
    - "LowerContinentalCrustStrongZone"
- Parameters:
    - `$(PARAMS.strain_final.name)` 
    - `$(PARAMS.friction_angle_initial.name)`
    - `$(PARAMS.friction_angle_final.name)`
## Radiogenic Heat Production in Continental Crust
- type_names:
    - "FelsicContinentalCrustFertile"
    - "FelsicContinentalCrustPartiallyMolten"
    - "FelsicContinentalCrustRefactory"
    - "MaficContinentalCrustFertile"
    - "MaficContinentalCrustPartiallyMolten"
    - "MaficContinentalCrustRefactory"
- Parameters:
    - `$(PARAMS.heat_production_upper_crust.name)`
    - `$(PARAMS.heat_production_lower_crust.name)`
## Pre-Exponential Factors in Asthenosphere
- type_names:
    If not in domains ["UpperMantleLithosphere", "MiddleMantleLithosphere", "LowerMantleLithosphere"]):
        - "UltramaficMantleFertile"
        - "UltramaficMantlePartiallyMolten"
        - "UltramaficMantleRefactory"
- Parameters:
    - `$(PARAMS.scale_factor_mantle_dislocation_creep.name)`
    - `$(PARAMS.scale_factor_mantle_diffusion_creep.name)`
## Pre-Exponential Factors in Continental Crust
- type_names:
    - "FelsicContinentalCrustFertile"
    - "FelsicContinentalCrustPartiallyMolten"
    - "FelsicContinentalCrustRefactory"
    - "MaficContinentalCrustFertile"
    - "MaficContinentalCrustPartiallyMolten"
    - "MaficContinentalCrustRefactory"
- Parameters:
    - `$(PARAMS.scale_factor_crustal_dislocation_creep.name)`
    - `$(PARAMS.scale_factor_oceanic_crust_dislocation_creep.name)`
## Pre-Exponential Factors in Oceanic Crust
- type_names:
    - "SolidifiedGabbro"
    - "SolidifiedGabbroPartiallyMolten"
    - "SolidifiedLayeredGabbro"
    - "SolidifiedLayeredGabbroPartiallyMolten"
    - "SolidifiedBasalt"
- Parameters:
    - `$(PARAMS.scale_factor_oceanic_crust_dislocation_creep.name)`
## Pre-Exponential Factors in Solidified Basalt
- type_names: "SolidifiedBasalt"
- Parameters:
    - `$(PARAMS.scale_factor_solidified_basalt_dislocation_creep.name)`
## Pre-Exponential Factors in Sediment
- type_names: "Sediment"
- Parameters:
    - `$(PARAMS.scale_factor_sediment_dislocation_creep.name)`
## Latent Heat in Mantle Rocks
- type_names:
    - "UltramaficMantleFertile"
    - "UltramaficMantlePartiallyMolten"
    - "UltramaficMantleRefactory"
- Parameters:
    - `$(PARAMS.latent_heat_mantle.name)`
## Latent Heat in Oceanic Crust
- type_names:
    - "SolidifiedGabbro"
    - "SolidifiedGabbroPartiallyMolten"
    - "SolidifiedLayeredGabbro"
    - "SolidifiedLayeredGabbroPartiallyMolten"
    - "SolidifiedBasalt"
- Parameters:
    - `$(PARAMS.latent_heat_oceanic_crust.name)`
## Mantle Solidus
- type_names:
    - "UltramaficMantleFertile"
    - "UltramaficMantlePartiallyMolten"
    - "UltramaficMantleRefactory"
- Parameters:
    - `$(PARAMS.mantle_solidus.name)`
## Mantle Liquidus
- type_names:
    - "UltramaficMantleFertile"
    - "UltramaficMantlePartiallyMolten"
    - "UltramaficMantleRefactory"
- Parameters:
    - $(PARAMS.mantle_liquidus.name)

"""
function override_material_properties!(model::ModelData, case_inputs_active::CaseType)
    input_dict = Dict(Symbol(key) => case_inputs_active[key].value for key in Base.keys(case_inputs_active))
    override_material_properties(model; input_dict)
end

""" Override material properties for model materials. 

Note that the override will only occur if the keyword argument is provided. If 
the keyword argument is not provided, the material property will not be overridden.
"""
function override_material_properties(model::ModelData; input_dict::Dict{Symbol, Any})
    override_plasticity_in_zones_with_strain_weakening(model; input_dict)
    override_plasticity_in_strong_zones(model; input_dict)
    override_plasticity_in_oceanic_crust(model; input_dict)
    override_plasticity_in_solidified_basalt(model; input_dict)
    override_plasticity_in_sediment(model; input_dict)
    override_pre_exponential_factors_in_asthenosphere(model; input_dict)
    override_pre_exponential_factors_in_continental_crust(model; input_dict)
    override_pre_exponential_factors_in_oceanic_crust(model; input_dict)
    override_pre_exponential_factors_in_solidified_basalt(model; input_dict)
    override_pre_exponential_factors_in_sediment(model; input_dict)
    override_radiogenic_heat_production_in_continental_crust(model; input_dict)
    override_latent_heat_in_mantle_rocks(model; input_dict)
    override_latent_heat_in_oceanic_crust(model; input_dict)
    override_mantle_solidus(model; input_dict)
    override_mantle_liquidus(model; input_dict)
end

function override_plasticity_in_solidified_basalt(model::ModelData; input_dict::Dict{Symbol, Any})
    strain_final = 
        get(input_dict, Symbol(keys.strain_final_solidified_basalt.name), nothing)
    friction_angle_initial = 
        get(input_dict, Symbol(keys.friction_angle_initial_solidified_basalt.name), nothing)
    friction_angle_final = 
        get(input_dict, Symbol(keys.friction_angle_final_solidified_basalt.name), nothing)
    material_ids_for_basalt = get_ids_for_solidified_basalt(model)
    MaterialContainer.override_plasticity!(
        model.materials, strain_final, 
        friction_angle_initial, friction_angle_final, material_ids_for_basalt
    )
end

function override_plasticity_in_oceanic_crust(model::ModelData; input_dict::Dict{Symbol, Any})
    strain_final = 
        get(input_dict, Symbol(keys.strain_final_oceanic_crust.name), nothing)
    friction_angle_initial = 
        get(input_dict, Symbol(keys.friction_angle_initial_oceanic_crust.name), nothing)
    friction_angle_final = 
        get(input_dict, Symbol(keys.friction_angle_final_oceanic_crust.name), nothing)
    material_ids_for_oceanic_crust = get_ids_for_oceanic_crust(model)
    MaterialContainer.override_plasticity!(
        model.materials, strain_final, 
        friction_angle_initial, friction_angle_final, 
        material_ids_for_oceanic_crust
    )
end

function override_plasticity_in_sediment(model::ModelData; input_dict::Dict{Symbol, Any})
    strain_final = 
        get(input_dict, Symbol(keys.strain_final_sediment.name), nothing)
    friction_angle_initial = 
        get(input_dict, Symbol(keys.friction_angle_initial_sediment.name), nothing)
    friction_angle_final = 
        get(input_dict, Symbol(keys.friction_angle_final_sediment.name), nothing)
    material_ids_for_sediment = get_ids_for_sediment(model)
    MaterialContainer.override_plasticity!(
        model.materials, strain_final, 
        friction_angle_initial, friction_angle_final, 
        material_ids_for_sediment
    )
end

function override_plasticity_in_zones_with_strain_weakening(model::ModelData; input_dict::Dict{Symbol, Any})
    strain_final = 
        get(input_dict, Symbol(keys.strain_final.name), nothing)
    friction_angle_initial = 
        get(input_dict, Symbol(keys.friction_angle_initial.name), nothing)
    friction_angle_final = 
        get(input_dict, Symbol(keys.friction_angle_final.name), nothing)
    material_ids_to_update_plasticity = get_ids_for_plastic_materials_with_strain_weakening(model)
    MaterialContainer.override_plasticity_in_zones_with_strain_weakening!(
        model.materials, strain_final, 
        friction_angle_initial, friction_angle_final, 
        material_ids_to_update_plasticity
    )
end

function override_plasticity_in_strong_zones(model::ModelData; input_dict::Dict{Symbol, Any})
    strain_final = 
        get(input_dict, Symbol(keys.strain_final.name), nothing)
    friction_angle_initial = 
        get(input_dict, Symbol(keys.friction_angle_initial.name), nothing)
    friction_angle_final = 
        get(input_dict, Symbol(keys.friction_angle_final.name), nothing)
    material_ids_strong_zones = get_ids_for_lithospheric_strong_zones(model)
    MaterialContainer.override_plasticity_in_strong_zones!(
        model.materials, strain_final, 
        friction_angle_initial, friction_angle_final, 
        material_ids_strong_zones
    )
end

function override_radiogenic_heat_production_in_continental_crust(
    model::ModelData; 
    input_dict::Dict{Symbol, Any}
)
    heat_production_upper_crust = 
        get(input_dict, Symbol(keys.heat_production_upper_crust.name), nothing)
    heat_production_lower_crust = 
        get(input_dict, Symbol(keys.heat_production_lower_crust.name), nothing)
    (
        material_ids_felsic_continental_crust, material_ids_mafic_continental_crust
    ) = get_ids_for_felsic_and_mafic_continental_crust(model)
    MaterialContainer.override_radiogenic_heat_production_in_continental_crust!(
        model.materials, 
        heat_production_upper_crust,
        heat_production_lower_crust,
        material_ids_felsic_continental_crust,
        material_ids_mafic_continental_crust
    )
end

function override_latent_heat_in_mantle_rocks(model::ModelData; input_dict::Dict{Symbol, Any})
    latent_heat_new = 
        get(input_dict, Symbol(keys.latent_heat_mantle.name), nothing)
    material_ids_target = get_ids_for_all_mantle_rocks(model)
    MaterialContainer.override_latent_heat!(
        model.materials, latent_heat_new, material_ids_target
    )
end

function override_latent_heat_in_oceanic_crust(model::ModelData; input_dict::Dict{Symbol, Any})
    latent_heat_new = 
        get(input_dict, Symbol(keys.latent_heat_oceanic_crust.name), nothing)
    material_ids_target = get_ids_for_oceanic_crust(model)
    MaterialContainer.override_latent_heat!(
        model.materials, latent_heat_new, material_ids_target
    )
end

function override_mantle_solidus(model::ModelData; input_dict::Dict{Symbol, Any})
    solidus_model_new = 
        get(input_dict, Symbol(keys.mantle_solidus.name), nothing)
    material_ids_target = get_ids_for_all_mantle_rocks(model)
    MaterialContainer.override_solidus!(
        model.materials, String(solidus_model_new), material_ids_target
    )
end

function override_mantle_liquidus(model::ModelData; input_dict::Dict{Symbol, Any})
    liquidus_model_new = 
        get(input_dict, Symbol(keys.mantle_liquidus.name), nothing)
    material_ids_target = get_ids_for_all_mantle_rocks(model)
    MaterialContainer.override_liquidus!(
        model.materials, String(liquidus_model_new), material_ids_target
    )
end

function override_pre_exponential_factors_in_asthenosphere(model::ModelData; input_dict::Dict{Symbol, Any})
    scale_factor_mantle_dislocation_creep = 
        get(input_dict, Symbol(keys.scale_factor_mantle_dislocation_creep.name), nothing)
    scale_factor_mantle_diffusion_creep = 
        get(input_dict, Symbol(keys.scale_factor_mantle_diffusion_creep.name), nothing)
    material_ids_asthenosphere = get_ids_for_asthenospheric_mantle(model)
    MaterialContainer.override_pre_exponential_factors_in_asthenosphere!(
        model.materials, scale_factor_mantle_dislocation_creep,
        scale_factor_mantle_diffusion_creep, material_ids_asthenosphere
    )
end

function override_pre_exponential_factors_in_continental_crust(model::ModelData; input_dict::Dict{Symbol, Any})
    scale_factor_crust_dislocation_creep = 
        get(input_dict, Symbol(keys.scale_factor_crustal_dislocation_creep.name), nothing)
    material_ids_felsic_continental_crust, material_ids_mafic_continental_crust = 
        get_ids_for_felsic_and_mafic_continental_crust(model)
    material_ids_continental_crust = vcat(material_ids_felsic_continental_crust, material_ids_mafic_continental_crust)
    MaterialContainer.override_pre_exponential_factors_in_continental_crust!(
        model.materials, scale_factor_crust_dislocation_creep,
        material_ids_continental_crust
    )
end

function override_pre_exponential_factors_in_oceanic_crust(model::ModelData; input_dict::Dict{Symbol, Any})
    scale_factor_for_dislocation_creep = 
        get(input_dict, Symbol(keys.scale_factor_oceanic_crust_dislocation_creep.name), nothing)
    material_ids_target = get_ids_for_oceanic_crust(model)
    MaterialContainer.override_pre_exponential_factors_dislocation_creep!(
        model.materials, scale_factor_for_dislocation_creep, material_ids_target
    )
end

function override_pre_exponential_factors_in_solidified_basalt(model::ModelData; input_dict::Dict{Symbol, Any})
    scale_factor_for_dislocation_creep = 
        get(input_dict, Symbol(keys.scale_factor_solidified_basalt_dislocation_creep.name), nothing)
    material_ids_target = get_ids_for_solidified_basalt(model)
    MaterialContainer.override_pre_exponential_factors_dislocation_creep!(
        model.materials, scale_factor_for_dislocation_creep, material_ids_target
    )
end

function override_pre_exponential_factors_in_sediment(model::ModelData; input_dict::Dict{Symbol, Any})
    scale_factor_for_dislocation_creep = 
        get(input_dict, Symbol(keys.scale_factor_sediment_dislocation_creep.name), nothing)
    material_ids_target = get_ids_for_sediment(model)
    MaterialContainer.override_pre_exponential_factors_dislocation_creep!(
        model.materials, scale_factor_for_dislocation_creep, material_ids_target
    )
end

end # module 