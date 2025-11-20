module ElasticForecast

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: check_in_domain
using ....YieldStress: calc_yield_stress
using ....CheckYield: check_yield_is_applicable
using ..YieldCheck: check_for_yielding
using ..Limits: set_yield_stress_limits, apply_viscosity_limits

""" Update marker arrays for plastic yielding (itype_plasticity=1).

# Arguments
- `model`: Main model container object
- `no_yielding_in_mobile_wall`: Boolean to prevent yielding in mobile wall

# Updated Arrays
- `model.markers.arrays.rheology.marker_eta`: Marker viscosity (Pa.s)
- `model.markers.arrays.rheology.marker_pfailure`: Flag indicating plastic failure of marker
- `model.markers.arrays.stress.marker_sxx`: Marker normal stress (Pa)
- `model.markers.arrays.stress.marker_sxy`: Marker shear stress (Pa)

# Notes
- Marker stress arrays are only updated if the viscoelastic plastic yielding
  option (itype_plasticity = 0) is selected.
"""
function update_plastic_yielding!(
    model::ModelData,
    inside_flags::Vector{Int8},
    no_yielding_in_mobile_wall::Bool
)::Nothing
    # Unpack data structures for better performance
    iuse_fluid_pressure_for_yield = model.materials.parameters.stress_limits_yield.iuse_fluid_pressure_for_yield.value
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value

    viscosity_minimum = model.materials.parameters.viscosity_limits.viscosity_min.value
    viscosity_maximum = model.materials.parameters.viscosity_limits.viscosity_max.value
    stress_yield_minimum = model.materials.parameters.stress_limits_yield.yield_stress_min.value
    stress_yield_maximum = model.materials.parameters.stress_limits_yield.yield_stress_max.value

    marker_y = model.markers.arrays.location.marker_y.array

    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_sxx = model.markers.arrays.stress.marker_sxx.array
    marker_sxy = model.markers.arrays.stress.marker_sxy.array
    marker_exx = model.markers.arrays.strain.marker_exx.array
    marker_exy = model.markers.arrays.strain.marker_exy.array
    marker_pr = model.markers.arrays.pressure.marker_pr.array
    marker_sr_ratio = model.markers.arrays.strain.marker_sr_ratio.array
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_cohesion = model.markers.arrays.rheology.marker_cohesion.array
    marker_fric = model.markers.arrays.rheology.marker_fric.array
    marker_pfailure = model.markers.arrays.rheology.marker_pfailure.array

    mat_mu = model.materials.arrays.mat_mu.array

    domains = model.materials.dicts.matid_domains
    matid_sticky_air = domains["Atmosphere"]
    matid_mobile_wall = domains["MobileWall"]
    matid_plate_extension = domains["PlateExtension"]

    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                y_marker = marker_y[imarker]
                mat_id = marker_matid[imarker]
                cohesion = marker_cohesion[imarker]
                sine_friction_angle = marker_fric[imarker]
            end
            yield_is_applicable = check_yield_is_applicable(ntimestep, cohesion, sine_friction_angle)
            yielding_flag = 0
            if yield_is_applicable
                @inbounds begin
                    # Start with viscosity previously calculated from flow law
                    viscosity = marker_eta[imarker]
                    shear_modulus = mat_mu[mat_id]
                    stress_xx = marker_sxx[imarker]
                    stress_xy = marker_sxy[imarker]
                    strain_rate_xx = marker_exx[imarker]
                    strain_rate_xy = marker_exy[imarker]
                    strain_rate_ratio = marker_sr_ratio[imarker]
                    pressure = marker_pr[imarker]
                end
                stress_yield = calc_yield_stress(
                    cohesion, sine_friction_angle, pressure,
                    iuse_fluid_pressure_for_yield, y_marker, y_sealevel
                )
                stress_yield = set_yield_stress_limits(
                    stress_yield_minimum, stress_yield_maximum, stress_yield)
                # Plasticity yielding model from Gerya (2019) using purely
                # elastic stress forecast for yielding
                stress_invariant_forecast = forecast_marker_stress_elastic(
                    timestep, shear_modulus, stress_xx,
                    stress_xy, strain_rate_xx, strain_rate_xy,
                    strain_rate_ratio
                )
                (
                    yielding,
                    yielding_flag
                ) = check_for_yielding(
                    no_yielding_in_mobile_wall,
                    mat_id,
                    matid_sticky_air,
                    matid_mobile_wall,
                    matid_plate_extension,
                    stress_yield,
                    stress_invariant_forecast
                )
                if yielding
                    viscosity = correct_viscosity_for_yielding_itype1(
                        timestep,
                        shear_modulus,
                        stress_invariant_forecast,
                        stress_yield
                    )
                    viscosity = apply_viscosity_limits(
                        viscosity_minimum, viscosity_maximum, viscosity)
                    @inbounds marker_eta[imarker] = viscosity
                end
            end
            @inbounds begin
                if marker_pfailure[imarker] != 1.0f0
                    marker_pfailure[imarker] = Float32(yielding_flag)
                end
            end
        end
    end
    return nothing
end

@inline function forecast_marker_stress_elastic(
    timestep::Float64,
    shear_modulus::Float64,
    stress_xx::Float64,
    stress_xy::Float64,
    strain_rate_xx::Float64,
    strain_rate_xy::Float64,
    strain_rate_ratio::Float64
)::Float64
    stress_xx_new = (
        stress_xx
        + 2.0*shear_modulus*timestep*strain_rate_xx*strain_rate_ratio
    )
    stress_xy_new = (
        stress_xy
        + 2.0*shear_modulus*timestep*strain_rate_xy*strain_rate_ratio
    )

    stress_invariant = sqrt(stress_xx_new^2.0 + stress_xy_new^2.0)

    return stress_invariant
end

@inline function correct_viscosity_for_yielding_itype1(
    timestep::Float64,
    shear_modulus::Float64,
    stress_invariant_forecast::Float64,
    stress_yield::Float64
)::Float64
    viscosity = (
        shear_modulus*timestep*stress_yield
        / (stress_invariant_forecast - stress_yield)
    )

    return viscosity
end

end # module