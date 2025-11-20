module HeatCapacityTest

using Plots
import EarthBox.RockProperties.RhoCpModel.UpdateManager.TemperatureDependentWaples: 
    calculate_heat_capacity_waples
import EarthBox.ConversionFuncs: celsius_to_kelvin, kelvin_to_celsius

function run_test()::Nothing
    heat_capacity_at_20celsius = 800.0  # J/kg/K

    temperature_surface_celsius = 0.0
    density = 3300.0  # kg/m^3
    gravity = 9.8  # m/s^2

    ymax = 160_000.0  # m
    dy = 1000.0  # m
    ny = floor(Int, ymax/dy) + 1

    temperature_gradient = 1330/(ymax/1000.0)  # C/km

    depth_array = Float64[]
    pressure_array = Float64[]
    temperature_array = Float64[]
    heat_capacity_array = Float64[]

    for i in 1:ny
        y = (i - 1) * dy
        push!(depth_array, y/1000.0)

        temperature_celsius = temperature_surface_celsius + temperature_gradient*y/1000.0

        temperature_kelvin = (
            celsius_to_kelvin(temperature_surface_celsius)
            + temperature_gradient*y/1000.0
        )
        push!(temperature_array, kelvin_to_celsius(temperature_kelvin))

        pressure_pa = density*gravity*y
        push!(pressure_array, pressure_pa/1e9)

        heat_capacity = calculate_heat_capacity_waples(
            temperature_celsius, heat_capacity_at_20celsius
        )

        push!(heat_capacity_array, heat_capacity)
    end

    for depth_index in eachindex(depth_array)
        println("$(depth_index): Depth: $(depth_array[depth_index]), Pressure: $(pressure_array[depth_index]), Temperature: $(temperature_array[depth_index]), Heat Capacity: $(heat_capacity_array[depth_index])")
    end


    make_pressure_temperature_plot(
        depth_array, pressure_array,
        temperature_array, "TemperatureDependentWaples"
    )

    make_heat_capacity_plot(
        depth_array, heat_capacity_array, "TemperatureDependentWaples"
    )
    
    return nothing
end

function make_pressure_temperature_plot(
    depth_array::Vector{Float64},
    pressure_array::Vector{Float64},
    temperature_array::Vector{Float64},
    model_name::String
)::Nothing
    # Create the main plot with pressure
    p1 = plot(
        pressure_array, depth_array,
        label="Pressure",
        color=:red,
        xlabel="Pressure (GPa)",
        ylabel="Depth (km)",
        xlims=(0, 4.0),
        xticks=0:0.5:4,
        ylims=(0, 160.0),
        yticks=0:10:160,
        yaxis=:flip,
        size=(800, 600),
        legend=:bottomleft
    )
    
    # Add temperature on secondary x-axis
    plot!(
        twiny(),
        temperature_array, depth_array,
        label="Temperature",
        color=:blue,
        linestyle=:dot,
        xlabel="Temperature (C)",
        xlims=(0.0, 1400.0),
        xticks=0:100:1400,
        ylims=(0, 160.0),
        yticks=0:10:160,
        yaxis=:flip,
        legend=:topright
    )

    savefig(p1, "pressure_temperature_plot_$(model_name).png")

    return nothing
end

function make_heat_capacity_plot(
    depth_array::Vector{Float64},
    heat_capacity_array::Vector{Float64},
    model_name::String
)::Nothing
    # Create heat capacity plot
    p2 = plot(
        heat_capacity_array, depth_array,
        color=:black,
        label="Heat Capacity",
        xlabel="Heat Capacity (J/K/kg)",
        ylabel="Depth (km)",
        xlims=(600, 1500.0),
        xticks=600:100:1500,
        yaxis=:flip,
        yticks=0:10:160,
        size=(800, 600),
        legend=:bottomleft
    )
    savefig(p2, "heat_capacity_plot_$(model_name).png")

    return nothing
end

end # module 