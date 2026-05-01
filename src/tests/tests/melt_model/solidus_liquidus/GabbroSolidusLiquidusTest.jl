module GabbroSolidusLiquidusTest

using CairoMakie
import EarthBox.MeltModel.MeltFraction: get_melting_model_parameters
import EarthBox.ConversionFuncs: kelvin_to_celsius

function run_test()::Nothing
    plot_gabbro_solidus_liquidus()
    return nothing
end

function plot_gabbro_solidus_liquidus()::Nothing
    itype_solidus = 5  # 5 (Gerya 2010) or 6 (Schubert 2013)
    itype_liquidus = 5

    pressure_pascals = 2800 * 9.81 * 3000.0

    temperature_liquidus_list = Float64[]
    temperature_solidus_list = Float64[]
    pressure_list = Float64[]
    depth_km_list = Float64[]

    dz = 1000  # meters
    total_depth = 12000.0  # meters
    nvals = floor(Int, total_depth / dz) + 1

    for i in 1:nvals
        pressure_pascals = 2800 * 9.81 * ((i - 1) * dz)
        push!(pressure_list, pressure_pascals / 1e9)
        push!(depth_km_list, (i - 1) * dz / 1000.0)
        
        temperature_liquidus, temperature_solidus = get_melting_model_parameters(
            pressure_pascals,
            itype_solidus,
            itype_liquidus
        )
        push!(temperature_liquidus_list, kelvin_to_celsius(temperature_liquidus))
        push!(temperature_solidus_list, kelvin_to_celsius(temperature_solidus))
    end

    use_pressure = false
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "Temperature (C)")

    if use_pressure
        lines!(ax, temperature_liquidus_list, pressure_list; label = "Liquidus (Gerya 2010)")
        lines!(ax, temperature_solidus_list, pressure_list; label = "Solidus (Gerya 2010)")
        ax.ylabel = "Pressure_GPa"
    else
        lines!(ax, temperature_liquidus_list, depth_km_list; label = "Liquidus (Gerya 2010)")
        lines!(ax, temperature_solidus_list, depth_km_list; label = "Solidus (Gerya 2010)")
        ax.ylabel = "Depth (km)"
        ylims!(ax, 0, 12)
    end

    xlims!(ax, 950, 1500)
    ax.yreversed = true
    axislegend(ax)
    save("gabbro_solidus_liquidus.png", fig)
    return nothing
end

end # module 