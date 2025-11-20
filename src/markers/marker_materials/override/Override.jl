module Override

import EarthBox.ParameterRegistry: get_eb_parameters
import ..MaterialsContainer: Materials
import ..MaterialsContainer.MaterialProperties: set_value!
import ..GetMaterialIDs: get_ids_for_plastic_materials_with_strain_weakening
import ..GetMaterialIDs: get_ids_for_lithospheric_strong_zones
import ..GetMaterialIDs: get_ids_for_continental_crust
import ..GetMaterialIDs: get_ids_for_felsic_and_mafic_continental_crust
import ..GetMaterialIDs: get_ids_for_asthenospheric_mantle

function override_material_properties!(
    materials::Materials; 
    input_dict::Dict{Symbol, Any}
)::Nothing
    """ Override material properties at the material level.

    Note this does not update material arrays used in the model. It only
    updates the material properties in the materials object. This is used
    to override the material properties in plotting function where a model data
    object is not available.

    See model_data/materials/MaterialContainer.jl for the implementation that 
    updates material properties in the model arrays.
    """
    override_plasticity_in_zones_with_strain_weakening!(materials; input_dict)
    override_plasticity_in_strong_zones!(materials; input_dict)
    override_radiogenic_heat_production_in_continental_crust!(materials; input_dict)
    override_pre_exponential_factors_in_asthenosphere!(materials; input_dict)
    override_pre_exponential_factors_in_continental_crust!(materials; input_dict)
    return nothing
end

function override_plasticity_in_zones_with_strain_weakening!(
    materials::Materials;
    input_dict::Dict{Symbol, Any}
)::Nothing
    """ Override plasticity in zones with strain weakening.
    """
    keys = get_eb_parameters()
    (
        material_ids_zones_with_strain_weakening
    ) = get_ids_for_plastic_materials_with_strain_weakening(materials)
    materials_dict = materials.materials
    
    friction_angle_initial = get(input_dict, keys.friction_angle_initial.name, nothing)
    friction_angle_final = get(input_dict, keys.friction_angle_final.name, nothing)
    strain_final = get(input_dict, keys.strain_final.name, nothing)
    
    for material_id in material_ids_zones_with_strain_weakening
        material = materials_dict[string(material_id)]
        flow_law = material.flow_law
        plas = flow_law.plasticity
        set_value!(
            plas, plas.friction_angle_initial.name, 
            friction_angle_initial
        )
        set_value!(
            plas, plas.friction_angle_final.name, 
            friction_angle_final
        )
        set_value!(
            plas, plas.strain_final.name, 
            strain_final
        )
    end
    return nothing
end

function override_plasticity_in_strong_zones!(
    materials::Materials;
    input_dict::Dict{Symbol, Any}
)::Nothing
    """ Override plasticity in strong zones (zones without weakening).

    For strong zones the friction angle is constant.
    """
    keys = get_eb_parameters()
    material_ids_strong_zones = get_ids_for_lithospheric_strong_zones(materials)
    materials_dict = materials.materials
    
    friction_angle_initial = get(input_dict, keys.friction_angle_initial.name, nothing)
    friction_angle_final = get(input_dict, keys.friction_angle_final.name, nothing)
    strain_final = get(input_dict, keys.strain_final.name, nothing)
    
    for material_id in material_ids_strong_zones
        material = materials_dict[string(material_id)]
        flow_law = material.flow_law
        plas = flow_law.plasticity
        set_value!(
            plas, plas.friction_angle_initial.name, 
            friction_angle_initial
        )
        set_value!(
            plas, plas.friction_angle_final.name, 
            friction_angle_final
        )
        set_value!(
            plas, plas.strain_final.name, 
            strain_final
        )
    end
    return nothing
end

function override_pre_exponential_factors_in_asthenosphere!(
    materials::Materials;
    input_dict::Dict{Symbol, Any}
)::Nothing
    """ Override pre-exponential factors in asthenosphere.
    """
    keys = get_eb_parameters()
    material_ids_asthenosphere = get_ids_for_asthenospheric_mantle(materials)
    materials_dict = materials.materials
    
    scale_factor_mantle_dislocation_creep = get(
        input_dict, keys.scale_factor_mantle_dislocation_creep.name, nothing)
    scale_factor_mantle_diffusion_creep = get(
        input_dict, keys.scale_factor_mantle_diffusion_creep.name, nothing)
    
    for material_id in material_ids_asthenosphere
        material = materials_dict[string(material_id)]
        flow_law = material.flow_law
        
        if scale_factor_mantle_dislocation_creep !== nothing
            disloc_creep = flow_law.dislocation_creep
            pre_exponential = disloc_creep.pre_exponential_dc.value
            set_value!(
                disloc_creep,
                disloc_creep.pre_exponential_dc.name,
                pre_exponential * scale_factor_mantle_dislocation_creep
            )
        end
        
        if scale_factor_mantle_diffusion_creep !== nothing
            diff_creep = flow_law.diffusion_creep
            pre_exponential = diff_creep.pre_exponential_difc.value
            set_value!(
                diff_creep,
                diff_creep.pre_exponential_difc.name,
                pre_exponential * scale_factor_mantle_diffusion_creep
            )
        end
    end
    return nothing
end

function override_pre_exponential_factors_in_continental_crust!(
    materials::Materials;
    input_dict::Dict{Symbol, Any}
)::Nothing
    """ Override pre-exponential factors for dislocation creep in cont. crust.
    """
    keys = get_eb_parameters()
    scale_factor_crustal_dislocation_creep = get(
        input_dict, keys.scale_factor_crustal_dislocation_creep.name, nothing)
    if scale_factor_crustal_dislocation_creep === nothing
        return nothing
    end
    
    material_ids_continental_crust = get_ids_for_continental_crust(materials)
    materials_dict = materials.materials
    
    for material_id in material_ids_continental_crust
        material = materials_dict[string(material_id)]
        flow_law = material.flow_law
        disloc_creep = flow_law.dislocation_creep
        pre_exponential = disloc_creep.pre_exponential_dc.value
        set_value!(
            disloc_creep,
            disloc_creep.pre_exponential_dc.name,
            pre_exponential * scale_factor_crustal_dislocation_creep
        )
    end
    return nothing
end

function override_radiogenic_heat_production_in_continental_crust!(
    materials::Materials;
    input_dict::Dict{Symbol, Any}
)::Nothing
    """ Override radiogenic heat production in continental crust.
    """
    keys = get_eb_parameters()
    (
        material_ids_felsic_continental_crust,
        material_ids_mafic_continental_crust
    ) = get_ids_for_felsic_and_mafic_continental_crust(materials)
    materials_dict = materials.materials
    
    radiogenic_heat_production_felsic_crust = get(
        input_dict, keys.radiogenic_heat_production_felsic_crust.name, nothing)
    radiogenic_heat_production_mafic_crust = get(
        input_dict, keys.radiogenic_heat_production_mafic_crust.name, nothing)
    
    for material_id in material_ids_felsic_continental_crust
        material = materials_dict[string(material_id)]
        rhp = material.radiogenic_heatproduction
        set_value!(
            rhp,
            rhp.radiogenic_heat_production.name,
            radiogenic_heat_production_felsic_crust
        )
    end
    
    for material_id in material_ids_mafic_continental_crust
        material = materials_dict[string(material_id)]
        rhp = material.radiogenic_heatproduction
        set_value!(
            rhp,
            rhp.radiogenic_heat_production.name,
            radiogenic_heat_production_mafic_crust
        )
    end
    return nothing
end

end # module

