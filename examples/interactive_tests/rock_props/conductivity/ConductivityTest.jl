module ConductivityTest

using CairoMakie
import EarthBox.RockProperties.ThermalConductivityModel.UpdateManager.
    Liao14: thermal_conductivity_liao14
import EarthBox.RockProperties.ThermalConductivityModel.UpdateManager.SekiguchiWaples: 
    thermal_conductivity_sekiguchi_waples
import EarthBox.ConversionFuncs: celsius_to_kelvin, kelvin_to_celsius

const inputs_dict = Dict(
        "liao14" => Dict(
            # Ultra-basic rocks from Clauser and Huenges (1995) Table 1.
            "thermal_conductivity_k0" => 0.73,
            "thermal_conductivity_a" => 1293.0
        ),
        # Peridotite (typical) from Hantschel and Kauerauf (2009) Table E.4.
        "sekiguchi_waples" => Dict(
            "thermal_conductivity_k0" => 4.0
        )
    )

function run_test()::Nothing
    temperature_surface_celsius = 0.0
    density = 3300.0  # kg/m^3
    gravity = 9.8  # m/s^2

    ymax = 100_000.0  # m
    dy = 1000.0  # m
    ny = floor(Int, ymax/dy) + 1

    temperature_gradient = 1330/(ymax/1000.0)  # C/km

    for (model_name, input_list) in inputs_dict
        println("Running test for $(model_name) model...")

        depth_array = Float64[]
        pressure_array = Float64[]
        temperature_array = Float64[]
        conductivity_array = Float64[]

        for i in 1:ny
            y = (i - 1) * dy
            push!(depth_array, y/1000.0)

            temperature_kelvin = (
                celsius_to_kelvin(temperature_surface_celsius)
                + temperature_gradient*y/1000.0
            )
            push!(temperature_array, kelvin_to_celsius(temperature_kelvin))

            pressure_pa = density*gravity*y
            push!(pressure_array, pressure_pa/1e9)

            if model_name == "liao14"
                thermal_conductivity_k0 = input_list["thermal_conductivity_k0"]
                thermal_conductivity_a = input_list["thermal_conductivity_a"]
                thermal_conductivity = thermal_conductivity_liao14(
                    thermal_conductivity_k0, thermal_conductivity_a,
                    pressure_pa, temperature_kelvin
                )
            elseif model_name == "sekiguchi_waples"
                thermal_conductivity_k0 = input_list["thermal_conductivity_k0"]
                thermal_conductivity = thermal_conductivity_sekiguchi_waples(
                    thermal_conductivity_k0, temperature_kelvin
                )
            end

            push!(conductivity_array, thermal_conductivity)
        end

        make_pressure_temperature_plots(
            depth_array, pressure_array, temperature_array, model_name)

        make_conductivity_plot(
            depth_array, conductivity_array, model_name)

    end
    

    return nothing
end

function make_pressure_temperature_plots(
    depth_array::Vector{Float64},
    pressure_array::Vector{Float64},
    temperature_array::Vector{Float64},
    model_name::String
)::Nothing
    fig = Figure(size = (800, 600))
    ax_p = Axis(
        fig[1, 1];
        xlabel = "Pressure (GPa)",
        ylabel = "Depth (km)",
        xticks = 0:0.5:4,
        yticks = 0:10:100,
        xaxisposition = :bottom,
    )
    ax_t = Axis(
        fig[1, 1];
        xlabel = "Temperature (C)",
        xticks = 0:100:1400,
        xaxisposition = :top,
    )
    hidespines!(ax_t)
    hideydecorations!(ax_t)

    lines!(ax_p, pressure_array, depth_array;
           color = :red, label = "Pressure")
    lines!(ax_t, temperature_array, depth_array;
           color = :blue, linestyle = :dot, label = "Temperature")

    xlims!(ax_p, 0, 4.0)
    xlims!(ax_t, 0.0, 1400.0)
    ylims!(ax_p, 0, 100.0)
    ylims!(ax_t, 0, 100.0)
    ax_p.yreversed = true
    ax_t.yreversed = true

    axislegend(ax_p; position = :lb)
    axislegend(ax_t; position = :rt)

    save("pressure_temperature_plot_$(model_name).png", fig)

    return nothing
end

function make_conductivity_plot(
    depth_array::Vector{Float64},
    conductivity_array::Vector{Float64},
    model_name::String
)::Nothing

    fig = Figure(size = (800, 600))
    ax = Axis(
        fig[1, 1];
        xlabel = "Thermal Conductivity (W/m/K)",
        ylabel = "Depth (km)",
        yticks = 0:10:100,
    )
    lines!(ax, conductivity_array, depth_array;
           color = :black, label = "Thermal Conductivity")
    xlims!(ax, 1.5, 5.0)
    ax.yreversed = true
    axislegend(ax; position = :rb)
    save("conductivity_plot_$(model_name).png", fig)

    return nothing
end

end # module 
if abspath(PROGRAM_FILE) == @__FILE__
    ConductivityTest.run_test()
end
