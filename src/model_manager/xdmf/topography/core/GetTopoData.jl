module GetTopoData

import EarthBox.ModelDataContainer: ModelData
import ...XdmfUtils: intstr

function get_topo_data(model::ModelData)::Dict{String, Any}
    noutput = model.timestep.parameters.output_steps.noutput.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_Myr = model.conversion.parameters.sec_per_Myr.value

    model_time = timesum/sec_per_Myr
    time_units = "Myr"

    topo_x_m = copy(model.topography.arrays.gridt.array[1,:])
    topo_y_m = -copy(model.topography.arrays.gridt.array[2,:])

    jld_markerfile = "topo_" * intstr(noutput) * ".jld"

    ntopo = size(topo_x_m, 1)
    topo_xy_km = zeros(Float64, ntopo, 2)
    update_topo_xy_km_array(ntopo, topo_x_m, topo_y_m, topo_xy_km)
    topo_ids = collect(1:ntopo)

    topo_data = Dict{String, Any}(
        "ntopo" => ntopo,
        "jld_markerfile" => jld_markerfile,
        "jld_dataname_xy" => "topo_xy_km",
        "jld_dataname_ids" => "topo_ids",
        "topo_id_array" => topo_ids,
        "topo_xy_km_array" => topo_xy_km,
        "time" => model_time,
        "time_units" => time_units,
        "noutput" => noutput
    )

    return topo_data
end

function update_topo_xy_km_array(
    ntopo::Int,
    topo_x_m::Vector{Float64},
    topo_y_m::Vector{Float64},
    topo_xy_km::Matrix{Float64}
)::Nothing
    for i in 1:ntopo
        topo_xy_km[i, 1] = topo_x_m[i]/1000.0
        topo_xy_km[i, 2] = topo_y_m[i]/1000.0
    end
    return nothing
end

end # module 