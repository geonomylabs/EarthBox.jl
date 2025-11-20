module LoadMaterialArrays

import EarthBox.ModelDataContainer: ModelData
import EarthBox: OptionTools
import EarthBox.Rheology: FlowViscosity
import EarthBox.SolidusLiquidus: LiquidusModels
import EarthBox.SolidusLiquidus: SolidusModels
import ..MaterialContainer: Material
import ..MaterialProperties: set_flow_type!
import ..MaterialProperties: set_itype_solidus!
import ..MaterialProperties: set_itype_liquidus!

MaterialCollectionDictType = Dict{Union{String, Int16}, Material}

function load_material_arrays(materials::MaterialCollectionDictType, model::ModelData)
    matids = collect(keys(materials))
    for matid in matids
        if typeof(matid) == String
            matid = parse(Int16, matid)
        end
        if matid < 1
            error("matid must be greater than 0")
        end
        load_density_array(materials, matid, model)
        load_flow_arrays(materials, matid, model)
        load_shear_modulus_array(materials, matid, model)
        load_plasticity_array(materials, matid, model)
        load_heat_capacity_array(materials, matid, model)
        load_thermal_conductivity_array(materials, matid, model)
        load_heat_production_array(materials, matid, model)
        load_melting_arrays(materials, matid, model)
        load_compaction_array(materials, matid, model)
    end
end

""" Load material model information into density material array.

# Updated Arrays
- model.materials.arrays.mat_rho.array: Array((nmats,4), Float64)
    1: Standard density (kg/m^3)
    2: Thermal expansion (1/K)
    3: Compressibility (1/Pa)
    4: Melt density (kg/m^3)
"""
function load_density_array(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    density = materials[matid].density
    model.materials.arrays.mat_rho.array[matid,1] = density.standard_density.value
    model.materials.arrays.mat_rho.array[matid,2] = density.thermal_expansion.value
    model.materials.arrays.mat_rho.array[matid,3] = density.compressibility.value
    model.materials.arrays.mat_rho.array[matid,4] = density.melt_density.value
    return nothing
end

""" Load material model information into flow material arrays.

Updated Arrays:
- pymodel.materials.arrays.mat_flow_type.array: Array((nmats), Int64)
    Flow law type for each material:
    -1 = viscosity as a function of x
    0 = isoviscous
    1, 2 = composite rheology (dislocation, diffusion and Peierls creep)
    3 = not used
    4 = Temperature dependent viscosity for Couette flow benchmark
    5 = Temperature dependent viscosity for convection in a box benchmark

- pymodel.materials.arrays.mat_flow.array: Array((nmats, 17), Float64)
    1: Viscosity for isoviscous models (Pa.s)
    2: Pre-exponential factor for power law (1/s/MPa^n)
    3: Power law stress exponent n
    4: Activation energy for power law (kJ/mol)
    5: Activation volume for power law (J/MPa/mol)
    6: Pre-exponential term for Peierls creep (s^-1*MPa^-2)
    7: Peierls creep stress exponent 1
    8: Peierls creep stress exponent 2
    9: Peierls stress (MPa)
    10: Pre-exponential terms for diffusion creep (1/s/MPa^n)
    11: Activation energy for diffusion creep (kJ/mol)
    12: Activation volume for diffusion creep (J/MPa/mol)
    13: Pre-exponential term for temperature-dependent viscosity (Pa.s)
    14: Activation energy for temperature-dependent viscosity (kJ/mol)
    15: Reference viscosity for flow law from Blankenbach (1989) (Pa.s)
    16: b term for flow law from Blankenbach (1989)
    17: c term for flow law from Blankenbach (1989)
"""
function load_flow_arrays(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    flow_law = materials[matid].flow_law
    option_name = flow_law.flow_stype.value
    if option_name != "None"
        options = FlowViscosity.get_options()
        flow_type_value = OptionTools.get_id(options, option_name)
        set_flow_type!(flow_law, flow_type_value)
    end

    model.materials.arrays.mat_flow_type.array[matid] = Int(flow_law.flow_type.value)

    model.materials.arrays.mat_flow.array[matid,1] = flow_law.viscosity_iso.value
    model.materials.arrays.mat_flow.array[matid,2] = flow_law.dislocation_creep.pre_exponential_dc.value
    model.materials.arrays.mat_flow.array[matid,3] = flow_law.dislocation_creep.stress_exponent_n_dc.value
    model.materials.arrays.mat_flow.array[matid,4] = flow_law.dislocation_creep.activation_energy_dc.value
    model.materials.arrays.mat_flow.array[matid,5] = flow_law.dislocation_creep.activation_volume_dc.value
    model.materials.arrays.mat_flow.array[matid,6] = flow_law.peierls_creep.pre_exponential_pei.value
    model.materials.arrays.mat_flow.array[matid,7] = flow_law.peierls_creep.stress_exponent_m1_pei.value
    model.materials.arrays.mat_flow.array[matid,8] = flow_law.peierls_creep.stress_exponent_m2_pei.value
    model.materials.arrays.mat_flow.array[matid,9] = flow_law.peierls_creep.peierls_stress.value
    model.materials.arrays.mat_flow.array[matid,10] = flow_law.diffusion_creep.pre_exponential_difc.value
    model.materials.arrays.mat_flow.array[matid,11] = flow_law.diffusion_creep.activation_energy_difc.value
    model.materials.arrays.mat_flow.array[matid,12] = flow_law.diffusion_creep.activation_volume_difc.value
    model.materials.arrays.mat_flow.array[matid,13] = flow_law.temperature_dependent_viscosity.pre_exponential_td.value
    model.materials.arrays.mat_flow.array[matid,14] = flow_law.temperature_dependent_viscosity.activation_energy_td.value
    model.materials.arrays.mat_flow.array[matid,15] = flow_law.blankenbach89_viscosity.viscosity_ref_blankenbach89.value
    model.materials.arrays.mat_flow.array[matid,16] = flow_law.blankenbach89_viscosity.b_term_blankenbach89.value
    model.materials.arrays.mat_flow.array[matid,17] = flow_law.blankenbach89_viscosity.c_term_blankenbach89.value
    return nothing
end

""" Load material model information into shear modulus material array.

Updated Array:
- pymodel.materials.arrays.mat_mu.array: Shear modulus (Pa)
"""
function load_shear_modulus_array(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    shear_modulus = materials[matid].shear_modulus
    model.materials.arrays.mat_mu.array[matid] = shear_modulus.shear_modulus.value
    return nothing
end

""" Load material model information into plasticity material array.

Updated Array:
- pymodel.materials.arrays.mat_plastic.array: Array((nmats,7), Float64)
    1: Cohesion (C0) at initial strain (Pa)
    2: Cohesion (C1) at final strain (Pa)
    3: Sine of friction angle at initial strain
    4: Sine of friction angle at final strain
    5: Initial strain
    6: Final strain
    7: Dilatation angle
"""
function load_plasticity_array(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    plasticity = materials[matid].plasticity
    model.materials.arrays.mat_plastic.array[matid,1] = plasticity.cohesion_initial.value
    model.materials.arrays.mat_plastic.array[matid,2] = plasticity.cohesion_final.value
    model.materials.arrays.mat_plastic.array[matid,3] = plasticity.sine_friction_angle_initial.value
    model.materials.arrays.mat_plastic.array[matid,4] = plasticity.sine_friction_angle_final.value
    model.materials.arrays.mat_plastic.array[matid,5] = plasticity.strain_initial.value
    model.materials.arrays.mat_plastic.array[matid,6] = plasticity.strain_final.value
    model.materials.arrays.mat_plastic.array[matid,7] = plasticity.dilatation_angle.value
    return nothing
end

""" Load material model information into heat capacity material array.

Updated Array:
- pymodel.materials.arrays.mat_cp.array: Array((nmats), Float64)
    Heat capacity (J/kg/K)
"""
function load_heat_capacity_array(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    heat_capacity = materials[matid].heat_capacity
    model.materials.arrays.mat_cp.array[matid] = heat_capacity.heat_capacity.value
    return nothing
end

""" Load material model information into conductivity material array.

Updated Array:
- pymodel.materials.arrays.mat_kt.array: Array((nmats,2), Float64)
    1: Reference thermal conductivity k0 (W/m/k)
    2: Thermal conductivity term a (W/m)
"""
function load_thermal_conductivity_array(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    thermal_cond = materials[matid].thermal_conductivity
    model.materials.arrays.mat_kt.array[matid,1] = thermal_cond.thermal_conductivity_ref.value
    model.materials.arrays.mat_kt.array[matid,2] = thermal_cond.thermal_conductivity_a.value
    return nothing
end

""" Load material model info into heat production material array.

Updated Array:
- pymodel.materials.arrays.mat_hr.array: Array((nmats), Float64)
    Radiogenic heat production (W/m^3)
"""
function load_heat_production_array(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    heatprod = materials[matid].radiogenic_heatproduction
    model.materials.arrays.mat_hr.array[matid] = heatprod.radiogenic_heat_production.value
    return nothing
end

""" Load material model info into melt material arrays.

Updated Arrays:
- pymodel.materials.arrays.mat_melting_itypes.array: Array((nmats,2), Int64)
    1: Solidus itype
    2: Liquidus itype

- pymodel.materials.arrays.mat_melting.array: Array((nmats,2), Float64)
    1: Latent heat (J/kg)
    2: Solidus temperature shift to approximate wet conditions (K)
"""
function load_melting_arrays(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    melting_parameters = materials[matid].melting_parameters

    # Override default values for solidus and liquidus itypes
    option_name = melting_parameters.stype_solidus.value
    if option_name != "None"
        solidus_options = SolidusModels.get_options()
        solidus_id = OptionTools.get_id(solidus_options, option_name)
        set_itype_solidus!(melting_parameters, solidus_id)
    end

    option_name = melting_parameters.stype_liquidus.value
    if option_name != "None"
        liquidus_options = LiquidusModels.get_options()
        liquidus_id = OptionTools.get_id(liquidus_options, option_name)
        set_itype_liquidus!(melting_parameters, liquidus_id)
    end

    model.materials.arrays.mat_melting_itypes.array[matid,1] = melting_parameters.itype_solidus.value
    model.materials.arrays.mat_melting_itypes.array[matid,2] = melting_parameters.itype_liquidus.value
    model.materials.arrays.mat_melting.array[matid,1] = melting_parameters.latent_heat.value
    return nothing
end

""" Load material model information into compaction material array.

Updated Array:
- pymodel.materials.arrays.mat_compaction.array: Array((nmats,2), Float64)
    1: Initial porosity at sediment-water interface (fraction)
    2: Porosity decay depth (m) used in Athy's law
"""
function load_compaction_array(
    materials::MaterialCollectionDictType,
    matid::Int16,
    model::ModelData
)::Nothing
    compaction = materials[matid].compaction
    model.materials.arrays.mat_compaction.array[matid,1] = compaction.porosity_initial.value
    model.materials.arrays.mat_compaction.array[matid,2] = compaction.porosity_decay_depth.value
    return nothing
end

end # module LoadMaterialArrays 