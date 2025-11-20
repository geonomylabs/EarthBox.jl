module PeridotiteSolidusLiquidusTest

using Plots
import EarthBox.MeltModel.MeltFraction: get_melting_model_parameters
import EarthBox.ConversionFuncs: kelvin_to_celsius

const melting_models = Dict(
        "Gerya2010" => Dict("itype_solidus" => 3, "itype_liquidus" => 3),
        "Katz2003" => Dict("itype_solidus" => 4, "itype_liquidus" => 4)
    )

function run_test()::Nothing
    plot_peridotite_solidus_liquidus(model_name="Gerya2010")
    return nothing
end

function plot_peridotite_solidus_liquidus(;model_name::String)::Nothing
    # Model name and itype for solidus and liquidus
    itype_solidus = melting_models[model_name]["itype_solidus"]
    itype_liquidus = melting_models[model_name]["itype_liquidus"]

    pressure_pascals = 2800 * 9.81 * 3000.0

    temperature_liquidus_list = Float64[]
    temperature_solidus_list = Float64[]
    pressure_list = Float64[]
    depth_km_list = Float64[]

    dz = 5000  # meters
    total_depth = 660_000.0  # meters
    nvals = floor(Int, total_depth / dz) + 1

    for i in 1:nvals
        pressure_pascals = 3300 * 9.81 * ((i - 1) * dz)
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

    use_pressure = true
    p = plot()
    plot!(p, legend=:bottomleft)
    if use_pressure
        plot!(p, temperature_liquidus_list, pressure_list, label="Liquidus ($model_name)")
        plot!(p, temperature_solidus_list, pressure_list, label="Solidus ($model_name)")
        ylabel!("Pressure_GPa")
    else
        plot!(p, temperature_liquidus_list, depth_km_list, label="Liquidus ($model_name)")
        plot!(p, temperature_solidus_list, depth_km_list, label="Solidus ($model_name)")
        ylabel!("Depth (km)")
        ylims!(0, 6)
    end
    xticks!(1000:200:2300)
    xlims!(1000, 2300)
    xlabel!("Temperature (C)")
    yaxis!(:flip)
    savefig(p, "peridotite_solidus_liquidus_$(model_name).png")
    return nothing
end

end # module 