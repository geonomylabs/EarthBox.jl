module CompositeYieldStress

import EarthBox.Rheology.CompositeCreep.PeierlsCreep: calculate_creep_stress_peierls
import EarthBox.Rheology.CompositeCreep: use_peierls_creep
import EarthBox.Rheology.CompositeCreep.Creep: calculate_creep_stress
import EarthBox.Rheology.PlasticFailure.YieldStress: calc_yield_stress as calculate_plastic_yield_stress
import EarthBox.Rheology.CompositeCreep.ExponentialTerm: calculate_activation_energy_term
import ..GetRheology: Rheology
import ..GetBasicLithosphereIDs: LithosphereMaterialIDs

function calculate_lithosphere_composite_yield_stress(
    strain_rate::Float64,
    lith_ids::LithosphereMaterialIDs,
    matid_gridy::Vector{Int16},
    upr_crust_rheology::Rheology,
    lwr_crust_rheology::Rheology,
    mantle_lithosphere_rheology::Rheology,
    asthenosphere_rheology::Rheology,
    temp_gridy::Vector{Float64},
    pressure_gridy::Vector{Float64},
    gridy::Vector{Float64};
    use_damaged_state::Bool = false,
    gas_constant::Float64 = 8.314, # J/mol/K
    iuse_fluid_pressure_for_yield::Int = 0
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    id_upr_crust = lith_ids.upper_continental_crust
    id_lwr_crust = lith_ids.lower_continental_crust
    id_mantle_lithosphere = lith_ids.mantle_lithosphere
    id_asthenosphere = lith_ids.asthenosphere

    nnodes = length(matid_gridy)
    yield_stress = zeros(Float64, nnodes)
    dislocation_stress = zeros(Float64, nnodes)
    diffusion_stress = zeros(Float64, nnodes)
    peierls_stress = zeros(Float64, nnodes)

    for i in 1:nnodes
        y_location = gridy[i]
        material_id = matid_gridy[i]

        rheology = if material_id == id_upr_crust
            upr_crust_rheology
        elseif material_id == id_lwr_crust
            lwr_crust_rheology
        elseif material_id == id_mantle_lithosphere
            mantle_lithosphere_rheology
        elseif material_id == id_asthenosphere
            asthenosphere_rheology
        else
            error("Material ID not found.")
        end

        disloc_creep = rheology.dislocation_creep
        diff_creep = rheology.diffusion_creep
        pei_creep = rheology.peierls_creep
        plastic_failure = rheology.plastic_failure

        temperature_kelvin = temp_gridy[i]
        pressure_pascals = pressure_gridy[i]

        if disloc_creep.pre_exponential > 0.0
            disloc_creep_stress, _ = calculate_creep_stress(
                pressure_pascals, temperature_kelvin, gas_constant,
                disloc_creep.pre_exponential, disloc_creep.activation_energy,
                disloc_creep.activation_volume,
                false, 0.0,
                strain_rate,
                stress_exponent=disloc_creep.stress_exponent
            )
        else
            disloc_creep_stress = 1e39
        end

        if diff_creep.pre_exponential > 0.0
            diff_creep_stress, _ = calculate_creep_stress(
                pressure_pascals, temperature_kelvin, gas_constant,
                diff_creep.pre_exponential, diff_creep.activation_energy,
                diff_creep.activation_volume,
                false, 0.0,
                strain_rate,
                stress_exponent=1.0
            )
        else
            diff_creep_stress = 1e39
        end

        if use_peierls_creep(pei_creep.pre_exponential, temperature_kelvin)
            pei_creep_stress = calculate_creep_stress_peierls(
                pressure_pascals, temperature_kelvin, gas_constant,
                pei_creep.pre_exponential, pei_creep.activation_energy,
                pei_creep.activation_volume,
                false, 0.0,
                pei_creep.peierls_stress,
                strain_rate,
                pei_creep.stress_exponent_m,
                pei_creep.stress_exponent_n
            )
        else
            pei_creep_stress = 1e39
        end

        if !use_damaged_state
            sine_of_friction_angle = sin(deg2rad(plastic_failure.friction_angle_initial))
            cohesion = plastic_failure.cohesion_initial
        else
            sine_of_friction_angle = sin(deg2rad(plastic_failure.friction_angle_final))
            cohesion = plastic_failure.cohesion_final
        end

        plastic_failure_stress = calculate_plastic_yield_stress(
            cohesion,
            sine_of_friction_angle,
            pressure_pascals,
            iuse_fluid_pressure_for_yield,
            y_location,
            0.0
        )

        _yield_stress = min(
            disloc_creep_stress, diff_creep_stress,
            pei_creep_stress, plastic_failure_stress
        )

        yield_stress[i] = _yield_stress
        dislocation_stress[i] = disloc_creep_stress
        diffusion_stress[i] = diff_creep_stress
        peierls_stress[i] = pei_creep_stress
    end

    return yield_stress, dislocation_stress, diffusion_stress, peierls_stress
end

end # module CompositeYieldStress 