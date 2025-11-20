module LatentHeatTest

using Plots
import EarthBox.MeltModel.MeltPropertiesOpt: calculate_delta_meltfrac_times_latent_heat_for_temperature
import EarthBox.MeltModel.MeltPropertiesOpt: calculate_delta_meltfrac_times_latent_heat_for_pressure

function run_test()::Nothing
    latent_heat = 400.0  # kJ/kg
    itype_solidus = 4  # Katz 2003
    itype_liquidus = 4  # Katz 2003
    density = 3300.0  # kg/m^3
    latent_heat = latent_heat * 1.0e3  # J/kg
    delta_temperature = 1.0  # Kelvins
    delta_pressure = 1000.0  # Pa

    (
        pressure_pascals_list,
        depth_km_list,
        temperature_kelvin_list
    ) = get_geotherm()

    delta_heat_capacity_list = Float64[]
    delta_expansivity_list = Float64[]
    for i in eachindex(pressure_pascals_list)
        pressure_pascals = pressure_pascals_list[i]
        temperature_kelvins = temperature_kelvin_list[i]
        
        delta_xmelt_times_latent_heat = calculate_delta_meltfrac_times_latent_heat_for_temperature(
            pressure_pascals,
            temperature_kelvins,
            delta_temperature,
            itype_solidus,
            itype_liquidus,
            latent_heat
        )
        delta_heat_capacity = delta_xmelt_times_latent_heat / (2.0 * delta_temperature)
        push!(delta_heat_capacity_list, delta_heat_capacity)

        delta_xmelt_times_latent_heat = calculate_delta_meltfrac_times_latent_heat_for_pressure(
            delta_pressure,
            pressure_pascals,
            temperature_kelvins,
            itype_solidus,
            itype_liquidus,
            latent_heat
        )
        delta_adiabatic = density * delta_xmelt_times_latent_heat / (2.0 * delta_pressure)
        delta_expansivity = delta_adiabatic / temperature_kelvins
        push!(delta_expansivity_list, delta_expansivity * 1e5)
    end

    make_plot(
        delta_heat_capacity_list,
        delta_expansivity_list,
        depth_km_list,
        temperature_kelvin_list,
        pressure_pascals_list
    )

    return nothing
end

function get_geotherm()::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    pressure_pascals_list = Float64[]
    depth_km_list = Float64[]
    temp_kelvin_list = Float64[]
    
    degrees_per_km = 20.0
    dz = 1000  # meters
    total_depth = 140_000.0  # meters
    nvals = floor(Int, total_depth/dz) + 1
    
    for i in 1:nvals
        depth_m = (i - 1) * dz
        pressure_pascals = 3330 * 9.81 * depth_m
        push!(pressure_pascals_list, pressure_pascals)
        push!(depth_km_list, (i - 1) * dz / 1000.0)
        temperature_kelvin = 1330 + 273.15
        push!(temp_kelvin_list, temperature_kelvin)
    end
    
    return pressure_pascals_list, depth_km_list, temp_kelvin_list
end

function make_plot(
    delta_heat_capacity_list::Vector{Float64},
    delta_expansivity_list::Vector{Float64},
    depth_km_list::Vector{Float64},
    temperature_kelvin_list::Vector{Float64},
    pressure_pascals_list::Vector{Float64}
)::Nothing
    
    # Plot delta heat capacity
    p1 = plot(
        delta_heat_capacity_list, depth_km_list,
        xlabel="Delta Heat Capacity (J/kg/K)",
        ylabel="Depth (km)",
        yaxis=:flip
    )
    savefig(p1, "delta_heat_capacity.png")

    # Plot delta expansivity
    p2 = plot(
        delta_expansivity_list, depth_km_list,
        xlabel="Delta Expansivity (1/K) * 1e5",
        ylabel="Depth (km)",
        yaxis=:flip
    )
    savefig(p2, "delta_expansivity.png")

    # Plot temperature
    temperature_celsius_list = temperature_kelvin_list .- 273.15
    p3 = plot(
        temperature_celsius_list, depth_km_list,
        xlabel="Temperature (C)",
        ylabel="Depth (km)",
        yaxis=:flip
    )
    savefig(p3, "temperature.png")

    # Plot pressure
    pressure_gpa_list = pressure_pascals_list ./ 1e9
    p4 = plot(
        pressure_gpa_list, depth_km_list,
        xlabel="Pressure (GPa)",
        ylabel="Depth (km)",
        yaxis=:flip
    )
    savefig(p4, "pressure.png")

    return nothing
end

end # module 