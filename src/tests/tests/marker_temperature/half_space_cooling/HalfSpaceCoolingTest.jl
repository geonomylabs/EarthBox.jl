module HalfSpaceCoolingTest

using Plots
import EarthBox.Markers.MarkerTemperature.InitManager.HalfSpaceCooling: calculate_temperature
import EarthBox.ConversionFuncs: kelvin_to_celsius, celsius_to_kelvin,
    cm_per_yr_to_meters_per_seconds, seconds_to_years

function run_test()::Nothing
    surface_temperature_kelvin = celsius_to_kelvin(0.0)
    asthenosphere_temperature_kelvin = celsius_to_kelvin(1300.0)
    distance_from_ridge_meters = 1_000_000.0
    half_spreading_rate_m_s = cm_per_yr_to_meters_per_seconds(5.0)
    thermal_diffusivity_m2_s = 1.0e-6

    age_seconds = distance_from_ridge_meters / half_spreading_rate_m_s
    age_yrs = seconds_to_years(age_seconds)
    age_myr = age_yrs / 1.0e6

    dy = 500.0
    ymax = 100_000.0
    ny = floor(Int, ymax/dy) + 1

    profile_y = zeros(ny)
    temperature_profile_celsius = zeros(ny)
    for i in 1:ny
        depth_below_surface_meters = (i - 1) * dy
        profile_y[i] = depth_below_surface_meters
        temperature_kelvin = calculate_temperature(
            depth_below_surface_meters,
            surface_temperature_kelvin,
            asthenosphere_temperature_kelvin,
            distance_from_ridge_meters,
            half_spreading_rate_m_s,
            thermal_diffusivity_m2_s
        )
        temperature_profile_celsius[i] = kelvin_to_celsius(temperature_kelvin)
    end

    plot_profile(profile_y, temperature_profile_celsius, age_myr)

    return nothing
end

function plot_profile(
    profile_y::Vector{Float64},
    temperature_profile_celsius::Vector{Float64},
    age::Float64
)::Nothing
    
    p = plot(
        temperature_profile_celsius, profile_y ./ 1.0e3,
        xlabel="T (C)",
        ylabel="Depth (km)", 
        yaxis=:flip,
        title="Half-space cooling model at $(round(age, digits=2)) Myr",
        legend=:bottomleft
    )
    savefig(p, "half_space_profile.png")
    
    return nothing
end

end # module 