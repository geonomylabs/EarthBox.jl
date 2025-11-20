module GetRheology

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import ..GetBasicLithosphereIDs: LithosphereMaterialIDs

""" Creep parameters structure.

# Fields
- `pre_exponential::Float64`: Pre-exponential factor in 1/s1/MPa^n
- `activation_energy::Float64`: Activation energy in kJ/mol
- `activation_volume::Float64`: Activation volume in J/MPa/mol
- `stress_exponent::Float64`: Stress exponent n
"""
Base.@kwdef struct Creep
    pre_exponential::Float64 = 0.0
    activation_energy::Float64 = 0.0
    activation_volume::Float64 = 0.0
    stress_exponent::Float64 = 0.0
end

""" Peierls creep parameters structure.

# Fields
- `pre_exponential::Float64`: Pre-exponential factor in 1/s1/MPa^2
- `activation_energy::Float64`: Activation energy in kJ/mol
- `activation_volume::Float64`: Activation volume in J/MPa/mol
- `stress_exponent_m::Float64`: Stress exponent m
- `stress_exponent_n::Float64`: Stress exponent n
- `peierls_stress::Float64`: Peierls stress in MPa
"""
Base.@kwdef struct PeierlsCreep
    pre_exponential::Float64 = 0.0
    activation_energy::Float64 = 0.0
    activation_volume::Float64 = 0.0
    stress_exponent_m::Float64 = 0.0
    stress_exponent_n::Float64 = 0.0
    peierls_stress::Float64 = 0.0
end

""" Plastic failure parameters structure.

# Fields
- `friction_angle_initial::Float64`: Initial friction angle in degrees
- `friction_angle_final::Float64`: Final friction angle in degrees
- `cohesion_initial::Float64`: Initial cohesion in Pa
- `cohesion_final::Float64`: Final cohesion in Pa
"""
Base.@kwdef struct PlasticFailure
    friction_angle_initial::Float64 = 0.0
    friction_angle_final::Float64 = 0.0
    cohesion_initial::Float64 = 0.0
    cohesion_final::Float64 = 0.0
end

""" Rheology parameters structure.

# Fields
- `dislocation_creep::Creep`: Dislocation creep parameters
- `diffusion_creep::Creep`: Diffusion creep parameters
- `peierls_creep::PeierlsCreep`: Peierls creep parameters
- `plastic_failure::PlasticFailure`: Plastic failure parameters
"""
struct Rheology
    dislocation_creep::Creep
    diffusion_creep::Creep
    peierls_creep::PeierlsCreep
    plastic_failure::PlasticFailure
end

""" Get basic lithosphere rheology.

# Arguments
- `materials::Materials`: Materials object
- `lith_ids::LithosphereMaterialIDs`: Lithosphere material IDs

# Returns
- `Tuple{Rheology, Rheology, Rheology, Rheology}`: Tuple of rheology parameters for 
  upper crust, lower crust, mantle lithosphere, and asthenosphere
"""
function get_basic_lithosphere_rheology(
    materials::Materials,
    lith_ids::LithosphereMaterialIDs
)::Tuple{Rheology, Rheology, Rheology, Rheology}
    upr_crust_rheology = get_rheology(lith_ids.upper_continental_crust, materials)
    lwr_crust_rheology = get_rheology(lith_ids.lower_continental_crust, materials)
    mantle_lithosphere_rheology = get_rheology(lith_ids.mantle_lithosphere, materials)
    asthenosphere_rheology = get_rheology(lith_ids.asthenosphere, materials)
    
    return (
        upr_crust_rheology, lwr_crust_rheology,
        mantle_lithosphere_rheology, asthenosphere_rheology
    )
end

""" Make rheology.

# Arguments
- `material_id::Int16`: Material ID
- `materials::Materials`: Materials object

# Returns
- `Rheology`: Rheology parameters structure
"""
function get_rheology(material_id::Int16, materials::Materials)::Rheology
    materials_dict = materials.materials
    material = materials_dict[material_id]

    flow_law = material.flow_law

    disloc_creep = flow_law.dislocation_creep
    dislocation_creep = Creep(
        pre_exponential=disloc_creep.pre_exponential_dc.value,
        activation_energy=disloc_creep.activation_energy_dc.value,
        activation_volume=disloc_creep.activation_volume_dc.value,
        stress_exponent=disloc_creep.stress_exponent_n_dc.value
    )

    diff_creep = flow_law.diffusion_creep
    diffusion_creep = Creep(
        pre_exponential=diff_creep.pre_exponential_difc.value,
        activation_energy=diff_creep.activation_energy_difc.value,
        activation_volume=diff_creep.activation_volume_difc.value,
        stress_exponent=1.0
    )

    pei_creep = flow_law.peierls_creep
    peierls_creep = PeierlsCreep(
        pre_exponential=pei_creep.pre_exponential_pei.value,
        activation_energy=disloc_creep.activation_energy_dc.value,
        activation_volume=disloc_creep.activation_volume_dc.value,
        stress_exponent_m=pei_creep.stress_exponent_m1_pei.value,
        stress_exponent_n=pei_creep.stress_exponent_m2_pei.value,
        peierls_stress=pei_creep.peierls_stress.value
    )

    plas = flow_law.plasticity
    plastic_failure = PlasticFailure(
        friction_angle_initial=plas.friction_angle_initial.value,
        friction_angle_final=plas.friction_angle_final.value,
        cohesion_initial=plas.cohesion_initial.value,
        cohesion_final=plas.cohesion_final.value
    )

    rheology = Rheology(
        dislocation_creep,
        diffusion_creep,
        peierls_creep,
        plastic_failure
    )

    return rheology
end

end # module
