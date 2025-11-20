module ViscoelasticForecast

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: check_in_domain
using EarthBox.ViscoelasticFactor: calculate_viscoelastic_factor
using ....YieldStress: calc_yield_stress
using ....CheckYield: check_yield_is_applicable
using ..YieldCheck: check_for_yielding
using ..Limits: set_yield_stress_limits, apply_viscosity_limits

function update_plastic_yielding!(
    model::ModelData,
    inside_flags::Vector{Int8},
    no_yielding_in_mobile_wall::Bool
)::Nothing
    # Unpack data structures for better performance
    iuse_fluid_pressure_for_yield = model.materials.parameters.stress_limits_yield.iuse_fluid_pressure_for_yield.value
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value

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

    domains = model.materials.dicts.matid_domains
    matid_sticky_air = domains["Atmosphere"]
    matid_mobile_wall = domains["MobileWall"]
    matid_plate_extension = domains["PlateExtension"]

    mat_mu = model.materials.arrays.mat_mu.array

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
                # Plasticity yielding model from Gerya (2010) using a
                # viscoelastic stress forecast for yielding
                stress_invariant_forecast = forecast_marker_stress_viscoelastic(
                    viscosity, timestep, shear_modulus, stress_xx,
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
                    (
                        stress_xx, stress_xy
                    ) = correct_stress_for_yielding(
                        stress_xx, stress_xy, stress_yield)
                    viscosity = correct_viscosity_for_yielding(
                        strain_rate_xx, strain_rate_xy,
                        strain_rate_ratio, stress_yield
                    )
                    viscosity = apply_viscosity_limits(
                        viscosity_minimum, viscosity_maximum, viscosity)
                    @inbounds begin
                        marker_sxx[imarker] = stress_xx
                        marker_sxy[imarker] = stress_xy
                        marker_eta[imarker] = viscosity
                    end
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

@inline function forecast_marker_stress_viscoelastic(
    viscosity::Float64,
    timestep::Float64,
    shear_modulus::Float64,
    stress_xx::Float64,
    stress_xy::Float64,
    strain_rate_xx::Float64,
    strain_rate_xy::Float64,
    strain_rate_ratio::Float64
)::Float64
    viscoelastic_factor = calculate_viscoelastic_factor(
        viscosity, shear_modulus, timestep)
    stress_xx_new = (
        stress_xx*viscoelastic_factor
        + 2.0*viscosity*strain_rate_xx*strain_rate_ratio
        * (1.0 - viscoelastic_factor)
    )
    stress_xy_new = (
        stress_xy*viscoelastic_factor
        + 2.0*viscosity*strain_rate_xy*strain_rate_ratio
        * (1.0 - viscoelastic_factor)
    )
    stress_invariant = sqrt(stress_xx_new^2.0 + stress_xy_new^2.0)
    return stress_invariant
end

@inline function correct_stress_for_yielding(
    stress_xx::Float64,
    stress_xy::Float64,
    stress_yield::Float64
)::Tuple{Float64, Float64}
    stress_invariant_old = sqrt(stress_xx^2.0 + stress_xy^2.0)
    # This conditional statement was added to avoid a zero in the denominator
    if stress_invariant_old != 0.0
        stress_xx = stress_xx*stress_yield/stress_invariant_old
        stress_xy = stress_xy*stress_yield/stress_invariant_old
    else
        stress_xx = stress_yield
        stress_xy = stress_yield
    end
    return stress_xx, stress_xy
end

@inline function correct_viscosity_for_yielding(
    strain_rate_xx::Float64,
    strain_rate_xy::Float64,
    strain_rate_ratio::Float64,
    stress_yield::Float64
)::Float64
    # Maxwell model based strain rate
    strain_rate_invariant_old = (
        strain_rate_ratio
        * sqrt(strain_rate_xx^2.0 + strain_rate_xy^2.0)
    )
    viscosity = stress_yield/2.0/strain_rate_invariant_old
    return viscosity
end

@inline function print_yield_check_info(
    imarker::Int,
    x_marker::Float64,
    y_marker::Float64,
    mat_id::Int,
    cohesion::Float64,
    sine_of_friction_angle::Float64,
    yield_is_applicable::Bool
)::Nothing
    if imarker == 120100
        println(
            "Yielding at marker ", imarker, " ", x_marker, " ", y_marker,
            " of material ", mat_id,
            ": Cohesion:", cohesion,
            ": sin(theta):", sine_of_friction_angle,
            ": Yielding applicable:", yield_is_applicable
        )
    end
    return nothing
end

function print_yielding_info(
    imarker::Int,
    pressure::Float64,
    stress_yield::Float64,
    viscosity::Float64,
    shear_modulus::Float64,
    stress_xx::Float64,
    stress_xy::Float64,
    strain_rate_xx::Float64,
    strain_rate_xy::Float64,
    strain_rate_ratio::Float64,
    stress_invariant_forecast::Float64,
    yielding::Bool,
    yielding_flag::Int
)::Nothing
    if imarker == 120100
        println(
            "imarker:", imarker,
            " Pressure:", pressure,
            " Stress yield:", stress_yield,
            " Viscosity:", viscosity,
            " Shear modulus:", shear_modulus,
            " Stress xx:", stress_xx,
            " Stress xy:", stress_xy,
            " Strain rate xx:", strain_rate_xx,
            " Strain rate xy:", strain_rate_xy,
            " Strain rate ratio:", strain_rate_ratio,
            " Stress invariant forecast:", stress_invariant_forecast,
            " Yielding:", yielding,
            " Yielding flag:", yielding_flag
        )
    end
    return nothing
end

function print_yield_update(
    imarker::Int,
    imarker_target::Int,
    stress_xx::Float64,
    stress_xy::Float64,
    viscosity::Float64
)::Nothing
    if imarker == imarker_target
        println(
            "imarker:", imarker,
            " Stress xx:", stress_xx,
            " Stress xy:", stress_xy,
            " Viscosity:", viscosity
        )
    end
    return nothing
end

function print_stress_yield(
    imarker::Int,
    imarker_target::Int,
    pressure::Float64,
    cohesion::Float64,
    sine_friction_angle::Float64,
    stress_yield::Float64
)::Nothing
    if imarker == imarker_target
        println(
            "imarker:", imarker,
            " Pressure:", pressure,
            " Cohesion:", cohesion,
            " Sine of friction angle:", sine_friction_angle,
            " Stress yield:", stress_yield
        )
    end
    return nothing
end

function print_yield_stress_post_update(
    imarker::Int,
    imarker_target::Int,
    stress_yield_minimum::Float64,
    stress_yield_maximum::Float64,
    stress_yield::Float64
)::Nothing
    if imarker == imarker_target
        println(
            "imarker:", imarker,
            " stress_yield_min:",
            stress_yield_minimum,
            ": stress_yield_max:",
            stress_yield_maximum,
            ": stress_yield:",
            stress_yield
        )
    end
    return nothing
end

end # module