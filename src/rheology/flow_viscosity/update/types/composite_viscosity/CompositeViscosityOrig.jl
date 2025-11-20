module CompositeViscosity

using Printf
import EarthBox.ModelDataContainer: ModelData
import EarthBox.GridFuncs: check_in_domain
import EarthBox.Rheology.CompositeCreep: calculate_composite_creep_viscosity_and_strain_rate
import EarthBox.Rheology.ViscousStrainSoftening: update_pre_exponential_for_viscous_strain_softening

""" Calculate visco-elastic composite viscosity for a single marker

This function calculates the effective flow viscosity and pre-exponential term 
for each marker based on a composite viscosity model with visco-elastic stress 
forecast. Bisection is used to determine the effective viscosity that corresponds
to the forecasted stress invariant.

# Updated Arrays

model.markers.arrays.rheology.
------------------------------
marker_eta_flow.array: Vector{Float64}(undef, marknum)
    - Marker effective flow viscosity (Pa.s).

marker_preexp.array: Vector{Float64}(undef, marknum)
    - Power-law pre-exponential term (1/s/MPa^n).

# Inputs

model: ModelData
    - Structref model data object.

# Background

.. code-block::

Transient visco-elastic deviatoric stress is forecasted using a first-order
finite difference approximation as defined by the following equation:

sij = 2*eta_vp*eii*(1-ve_fac) + sij0*ve_fac                         (Eq. 1)

where sij is a deviatoric stress component, sij0 is the deviatoric stress
from the previous time step, eii is deviatoric strain rate, eta_vp is a
viscosity-like parameter equal to effective viscosity when no plastic
yielding occurs and ve_fac is a visco-elastic factor defined with the
following equation:

ve_fac = eta_vp / (mu*timestep + eta_vp)                            (Eq. 2)

Updated Stress and viscosity are calculated by applying bisection to
stress dependent flow laws and Equation 1 (see eq. 13.2 on page 179 of
Gerya, 2010).

Calculating stress-dependent effective viscosity when elastic deformation
is included requires an iterative approach since Equation (1) depends on
viscosity, which in turn depends on stress. Effective viscosity is
calculated using bisection where bisection limits are iteratively updated
as described in Figures (1) and (2).

                            Before Bisection Limit Update
                    |
                    |           x  sii_forecast
                    |
                    |         * sii_forecast_update
                    |
                sii |       o sii_avg = (sii_initial + sii_forecast)/2
                    |
                    |
                    |
                    |  x sii_initial (marker stress invariant)
                    |_______________
                            time

                            After Bisection Limit Update
                    |
                    |           x  sii_forecast (initial forecast)
                    |
                    |         * sii_avg (updated)
                    |
                sii |       x sii_initial (updated)
                    |
                    |
                    |
                    |
                    |_______________
                            time

    Figure 1: Bisection limit update for a case where stress forecast and
    updated stress forecast are greater than sii_avg.

The updated bisection limits for the case in Figure 1 are as follows:

    sii_initial = sii_avg

If the updated forecasted stress sii_forecast_updated is less than
sii_avg than bisection limits are updated as follows:

    sii_forecast = sii_avg


                            Before Bisection Limit Update
                    |
                    |           x  sii_forecast (initial forecast)
                    |
                    |
                    |
                sii |       o sii_avg = (sii_initial + sii_forecast)/2
                    |
                    |      * sii_forecast_update
                    |
                    |  x sii_initial (marker stress invariant)
                    |_______________
                            time

                            After Bisection Limit Update
                    |
                    |
                    |
                    |
                    |
                sii |       x sii_forecast (updated)
                    |
                    |     o sii_avg (updated)
                    |
                    |  x sii_initial (marker stress invariant)
                    |_______________
                            time

    Figure 2: Bisection limit update for a case where the updated stress
    forecast is less than sii_avg.

"""
function update_marker_flow_viscosity!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    timestep::Float64               = model.timestep.parameters.main_time_loop.timestep.value
    temperature_top::Float64        = model.bcs.parameters.temperature.temperature_top.value
    viscosity_minimum::Float64      = model.materials.parameters.viscosity_limits.viscosity_min.value
    viscosity_maximum::Float64      = model.materials.parameters.viscosity_limits.viscosity_max.value
    stress_minimum::Float64         = model.materials.parameters.stress_limits_power_law.powerlaw_stress_min.value
    gas_constant::Float64           = model.RGAS.value
    iuse_viscous_strain_soft::Int64 = model.materials.parameters.softening.iuse_viscous_strain_soft.value
    vsoftfac::Float64               = model.materials.parameters.softening.vsoftfac.value

    marknum::Int64 = model.markers.parameters.distribution.marknum.value

    marker_matid::Vector{Int16}      = model.markers.arrays.material.marker_matid.array
    marker_sxx::Vector{Float64}      = model.markers.arrays.stress.marker_sxx.array
    marker_sxy::Vector{Float64}      = model.markers.arrays.stress.marker_sxy.array
    marker_exx::Vector{Float64}      = model.markers.arrays.strain.marker_exx.array
    marker_exy::Vector{Float64}      = model.markers.arrays.strain.marker_exy.array
    marker_TK::Vector{Float64}       = model.markers.arrays.thermal.marker_TK.array
    marker_pr::Vector{Float64}       = model.markers.arrays.pressure.marker_pr.array
    marker_GII::Vector{Float64}      = model.markers.arrays.strain.marker_GII.array
    marker_sr_ratio::Vector{Float64} = model.markers.arrays.strain.marker_sr_ratio.array
    marker_eta_flow::Vector{Float64} = model.markers.arrays.rheology.marker_eta_flow.array
    marker_preexp::Vector{Float64}   = model.markers.arrays.rheology.marker_preexp.array

    mat_flow_type::Vector{Int64} = model.materials.arrays.mat_flow_type.array
    mat_flow::Matrix{Float64}    = model.materials.arrays.mat_flow.array
    mat_plastic::Matrix{Float64} = model.materials.arrays.mat_plastic.array
    mat_mu::Vector{Float64}      = model.materials.arrays.mat_mu.array

    max_iterations::Int64 = 20 # For bisection

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                mat_id = marker_matid[imarker]
                flow_type = mat_flow_type[mat_id]
            end
            if flow_type == 0
                @inbounds begin
                    viscosity = mat_flow[mat_id, 1]
                    pre_exponential_dislocation = mat_flow[mat_id, 2]
                end
            elseif flow_type == 1 || flow_type == 2
                # Marker information
                @inbounds begin
                    stress_xx = marker_sxx[imarker]
                    stress_xy = marker_sxy[imarker]
                    strain_rate_xx = marker_exx[imarker]
                    strain_rate_xy = marker_exy[imarker]
                    strain = marker_GII[imarker]
                    pressure = marker_pr[imarker]
                    strain_rate_ratio = marker_sr_ratio[imarker]
                    temperature = apply_minimum_temperature_limit(marker_TK[imarker], temperature_top)

                    # Material and flow law information
                    strain_initial = mat_plastic[mat_id, 5]
                    strain_final = mat_plastic[mat_id, 6]
                    pre_exponential_dislocation = mat_flow[mat_id, 2]
                    activation_energy_dislocation = mat_flow[mat_id, 4]
                    activation_volume_dislocation = mat_flow[mat_id, 5]
                    stress_exponent_dislocation = mat_flow[mat_id, 3]
                    pre_exponential_peierls = mat_flow[mat_id, 6]
                    stress_exponent_peierls_m = mat_flow[mat_id, 7]
                    stress_exponent_peierls_n = mat_flow[mat_id, 8]
                    stress_peierls = mat_flow[mat_id, 9]
                    pre_exponential_diffusion = mat_flow[mat_id, 10]
                    activation_energy_diffusion = mat_flow[mat_id, 11]
                    activation_volume_diffusion = mat_flow[mat_id, 12]
                    shear_modulus = mat_mu[mat_id]
                end

                pre_exponential_dislocation = apply_strain_softening(
                    iuse_viscous_strain_soft,
                    vsoftfac,
                    pre_exponential_dislocation,
                    strain,
                    strain_initial,
                    strain_final
                )

                stress_invariant_initial = calc_initial_stress_invariant(
                    stress_xx,
                    stress_xy,
                    stress_minimum
                )

                # Calculate initial effective viscosity using initial stress
                # invariant and static flow law terms so they can be reused
                (
                    viscosity,
                    _strain_rate,
                    exponential_term_dislocation,
                    activation_energy_term_peierls,
                    exponential_term_diffusion
                ) = calculate_composite_creep_viscosity_and_strain_rate(
                    pressure,
                    temperature,
                    gas_constant,
                    pre_exponential_dislocation,
                    activation_energy_dislocation,
                    activation_volume_dislocation,
                    stress_exponent_dislocation,
                    pre_exponential_peierls,
                    stress_exponent_peierls_m,
                    stress_exponent_peierls_n,
                    stress_peierls,
                    pre_exponential_diffusion,
                    activation_energy_diffusion,
                    activation_volume_diffusion,
                    stress_invariant_initial,
                    0, # iuse_old_term
                    0.0, # exponential_term_dislocation
                    0.0, # activation_energy_term_peierls
                    0.0, # exponential_term_diffusion
                )

                # Forecast second invariant of future marker stress using
                # initial viscosity
                stress_invariant_forecast = forecast_viscoelastic_stress(
                    timestep,
                    shear_modulus,
                    strain_rate_ratio,
                    stress_xx,
                    stress_xy,
                    strain_rate_xx,
                    strain_rate_xy,
                    stress_minimum,
                    viscosity
                )

                # Determine viscosity that corresponds to forecasted stress
                # invariant using bisection
                icount = 0
                delta_stress = 1e39
                while icount < max_iterations && delta_stress > 1.0
                    icount += 1
                    stress_invariant_avg = (
                        stress_invariant_initial + stress_invariant_forecast) / 2.0
                    (
                        viscosity, _, _, _, _
                    ) = calculate_composite_creep_viscosity_and_strain_rate(
                        pressure,
                        temperature,
                        gas_constant,
                        pre_exponential_dislocation,
                        activation_energy_dislocation,
                        activation_volume_dislocation,
                        stress_exponent_dislocation,
                        pre_exponential_peierls,
                        stress_exponent_peierls_m,
                        stress_exponent_peierls_n,
                        stress_peierls,
                        pre_exponential_diffusion,
                        activation_energy_diffusion,
                        activation_volume_diffusion,
                        stress_invariant_avg,
                        1, # iuse_old_term
                        exponential_term_dislocation,
                        activation_energy_term_peierls,
                        exponential_term_diffusion
                    )

                    stress_invariant_forecast_update = forecast_viscoelastic_stress(
                        timestep,
                        shear_modulus,
                        strain_rate_ratio,
                        stress_xx,
                        stress_xy,
                        strain_rate_xx,
                        strain_rate_xy,
                        stress_minimum,
                        viscosity
                    )

                    (
                        stress_invariant_initial,
                        stress_invariant_forecast
                    ) = update_bisection_limits(
                        stress_invariant_initial,
                        stress_invariant_forecast,
                        stress_invariant_avg,
                        stress_invariant_forecast_update
                    )

                    delta_stress = abs(stress_invariant_forecast - stress_invariant_initial)
                end

                viscosity = apply_viscosity_limits(
                    viscosity_minimum, viscosity_maximum, viscosity)
            end
            @inbounds begin
                marker_eta_flow[imarker] = viscosity
                marker_preexp[imarker] = pre_exponential_dislocation
            end
        end
    end

    return nothing
end

@inline function apply_strain_softening(
    iuse_viscous_strain_soft::Int,
    vsoftfac::Float64,
    pre_exponential_dislocation::Float64,
    strain::Float64,
    strain_initial::Float64,
    strain_final::Float64
)::Float64
    if iuse_viscous_strain_soft == 1
        pre_exponential_dislocation = update_pre_exponential_for_viscous_strain_softening(
            pre_exponential_dislocation,
            strain,
            vsoftfac,
            strain_initial,
            strain_final
        )
    end
    return pre_exponential_dislocation
end

@inline function calc_initial_stress_invariant(
    stress_xx::Float64,
    stress_xy::Float64,
    stress_minimum::Float64
)::Float64
    stress_invariant_initial = second_stress_invariant2d(stress_xx, stress_xy)
    return max(stress_invariant_initial, stress_minimum)
end

""" Forecast visco-elastic stress.

First-order finite difference yields the following equation for the
deviatoric stress transient:

sij = 2*eta_vp*eii*sr_ratio*(1-ve_fac) + sij0*ve_fac                (Eq. 1)

where sij is a deviatoric stress component, sij0 is the deviatoric stress
from the previous time step, eii is deviatoric strain rate, sr_ratio is
the ratio of strain rate calculated using grid stress changes and a
Maxwell model over strain rate interpolated from the grid, eta_vp is a
viscosity-like parameter equal to effective viscosity when no plastic
yielding occurs and ve_fac is a visco-elastic factor described with the
following equation:

    ve_fac = eta_vp / (mu*timestep + eta_vp)                        (Eq. 2)

Note that

    ve_fac = 1 - Z                                                  (Eq. 3)

where Z is defined in equation 13.2 on page 179 of Gerya (2010).

The strain-rate ratio sr_ratio is calculated by dividing the marker
strain rate calculated with grid stress changes and a Maxwell model by the
strain rate calculated with direct grid interpolation:

    sr_ratio = eii_marker_maxwell / eii_grid_interp_to_marker       (Eq. 4)

The second invariant of marker deviatoric strain rate eii in (Eq. 1) is
multiplied by sr_ratio to obtain a strain rate based on grid stress
changes:

    eii_dstress = eii*sr_ratio                                      (Eq. 5)

This approach reduces numerical diffusion. See pgs 190 and 191 of Gerya
(2010) for a more detailed explanation of why the strain-rate ratio is
used in (Eq. 1).

# Parameters

timestep : Float64
    Viscoelastic timestep (seconds)

shear_modulus : Float64
    Shear modulus (Pa) of marker material.

Viscosity : Float64
    Effective viscosity (Pa.s) of marker.

strain_rate_ratio : Float64
    Ratio of strain rate calculated using grid stress changes and a
    Maxwell model over strain rate interpolated from the grid.

stress_xx : Float64
    Deviatoric normal stress SIGMAxx (Pa).

stress_xy : Float64
    Deviatoric shear stress SIGMAxy (Pa).

strain_rate_xx : Float64
    Deviatoric normal strain rate EPSILONxx (1/s).

strain_rate_xy : Float64
    Deviatoric shear strain rate EPSILONxy (1/s).

# Returns
stress_invariant : Float64
    Forecasted second invariant of the deviatoric stress tensor.

"""
@inline function forecast_viscoelastic_stress(
    timestep::Float64,
    shear_modulus::Float64,
    strain_rate_ratio::Float64,
    stress_xx::Float64,
    stress_xy::Float64,
    strain_rate_xx::Float64,
    strain_rate_xy::Float64,
    stress_minimum::Float64,
    viscosity::Float64
)::Float64
    visco_elastic_factor = viscosity / (shear_modulus * timestep + viscosity)

    stress_xx_new = (
        2.0 * viscosity * strain_rate_xx * strain_rate_ratio * (1.0 - visco_elastic_factor) 
        + stress_xx * visco_elastic_factor
    )

    stress_xy_new = (
        2.0 * viscosity * strain_rate_xy * strain_rate_ratio * (1.0 - visco_elastic_factor) 
        + stress_xy * visco_elastic_factor
    )

    stress_invariant = second_stress_invariant2d(stress_xx_new, stress_xy_new)
    return max(stress_invariant, stress_minimum)
end

@inline function second_stress_invariant2d(stress_xx::Float64, stress_xy::Float64)::Float64
    return sqrt(stress_xx^2 + stress_xy^2)
end

@inline function update_bisection_limits(
    stress_invariant_initial::Float64,
    stress_invariant_forecast::Float64,
    stress_invariant_avg::Float64,
    stress_invariant_forecast_update::Float64
)::Tuple{Float64, Float64}
    if (stress_invariant_forecast > stress_invariant_initial && 
        stress_invariant_forecast_update > stress_invariant_avg) ||
       (stress_invariant_forecast < stress_invariant_initial && 
        stress_invariant_forecast_update < stress_invariant_avg)
        # Keep moving toward stress_invariant_forecast
        stress_invariant_initial = stress_invariant_avg
    else
        # Move toward stress_invariant_avg
        stress_invariant_forecast = stress_invariant_avg
    end
    return stress_invariant_initial, stress_invariant_forecast
end

@inline function apply_viscosity_limits(
    viscosity_minimum::Float64,
    viscosity_maximum::Float64,
    viscosity::Float64
)::Float64
    if viscosity < viscosity_minimum
        viscosity = viscosity_minimum
    elseif viscosity > viscosity_maximum
        viscosity = viscosity_maximum
    end
    return viscosity
end

@inline function apply_minimum_temperature_limit(
    temperature::Float64,
    temperature_minimum::Float64
)::Float64
    return max(temperature, temperature_minimum)
end

end # module 