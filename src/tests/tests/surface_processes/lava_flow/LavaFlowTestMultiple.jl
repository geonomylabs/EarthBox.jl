module LavaFlowTestMultiple

import Plots
import EarthBox.ConversionFuncs: mm_per_yr_to_meters_per_seconds as mm_yr_to_m_s
import EarthBox.ConversionFuncs: years_to_seconds
import EarthBox.ConversionFuncs: meters_per_year_to_meters_per_seconds as m_yr_to_m_s
import EarthBox.ConversionFuncs: meters_squared_per_year_to_meters_squared_per_second as m2_yr_to_m2_s
import EarthBox.ConversionFuncs: seconds_to_years
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.SurfaceProcesses.LavaFlowManager: calculate_number_of_flows
import EarthBox.SurfaceProcesses.LavaFlowManager: calculate_characteristic_volume_per_flow
import EarthBox.SurfaceProcesses.LavaFlowManager: forecast_eruption_location
import EarthBox.SurfaceProcesses.LavaFlowManager.LavaFlowSolverManager: LavaFlowSolver
import EarthBox.SurfaceProcesses.LavaFlowManager.LavaFlowSolverManager.MakeFlow.LavaFlowPulse: radiate_indices
import EarthBox.SurfaceProcesses.LavaFlowManager.LavaFlowSolverManager: extrude_magma
import EarthBox.TestManager.SedimentTransportTest: TopoGeometry, TopoGrids

function run_test()
    sediment_transport_parameters = SedimentTransportParameters(
        m2_yr_to_m2_s(0.25),
        m_yr_to_m_s(1.0), 
        1e-3,
        m2_yr_to_m2_s(1e2),
        2000.0,
        years_to_seconds(50_000.0),
        years_to_seconds(5.0*1e6),
        0.4,
        1/2000.0
    )

    topo_geometry = TopoGeometry(
        xsize=500_000.0, dx=500.0, xo_topo=50_000.0, width_topo=50_000.0,
        height_topo=-500.0, width_basin=300_000.0, depth_basin=0.0,
        y_sealevel=0.0, sediment_thickness_original=5000.0
    )

    decimation_factor = 10
    total_extrusion_volume = 1e7 # m3
    eruption_location_x_min = 200_000.0
    width_eruption_domain = 100_000.0
    residual_lava_thickness_subaerial = 20.0
    residual_laval_thickness_submarine = 50.0
    characteristic_flow_length_subaerial = 100_000.0
    characteristic_flow_length_submarine = 10_000.0

    (
        topo_gridx, topo_gridy_initial, _gridx_b, sediment_thickness_initial
    ) = initialize_grids(topo_geometry)

    eruption_location_x_forecast, eruption_location_y_forecast = 
        forecast_eruption_location(
            eruption_location_x_min, width_eruption_domain,
            topo_gridx, topo_gridy_initial
        )

    println("eruption_location_x_forecast: ", eruption_location_x_forecast)
    println("eruption_location_y_forecast: ", eruption_location_y_forecast)
    println("topo_geometry.y_sealevel: ", topo_geometry.y_sealevel)
    if eruption_location_y_forecast <= topo_geometry.y_sealevel
        println(">> eruption is forecasted to be subaerial")
    else
        println(">> eruption is forecasted to be submarine")
    end

    characteristic_volume_per_flow = calculate_characteristic_volume_per_flow(
        topo_geometry.y_sealevel, eruption_location_y_forecast,
        characteristic_flow_length_subaerial,
        characteristic_flow_length_submarine,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine
    )

    number_of_flows = calculate_number_of_flows(
        total_extrusion_volume, characteristic_volume_per_flow
    )

    println("   >> characteristic_volume_per_flow: ", characteristic_volume_per_flow)
    println("   >> number of flows: ", number_of_flows)

    lava_flow_solver = LavaFlowSolver(
        topo_gridx,
        topo_gridy_initial,
        sediment_thickness_initial,
        eruption_location_x_min,
        width_eruption_domain,
        total_extrusion_volume,
        number_of_flows,
        residual_lava_thickness_subaerial,
        residual_laval_thickness_submarine,
        topo_geometry.y_sealevel,
        sediment_transport_parameters,
        use_random_eruption_location=true,
        use_compaction_correction=true,
        decimation_factor=decimation_factor
    )

    extrude_magma(lava_flow_solver)

    top_sediment = lava_flow_solver.topo_gridy .+ lava_flow_solver.total_lava_thickness
    basement_initial = lava_flow_solver.topo_gridy_initial .+ sediment_thickness_initial

    make_plot(
        lava_flow_solver.topo_gridx,
        lava_flow_solver.topo_gridy,
        top_sediment,
        basement_initial
    )
end

function initialize_grids(
    topo_geometry::TopoGeometry
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    xsize = topo_geometry.xsize
    dx = topo_geometry.dx
    xnum = floor(Int, xsize/dx) + 1
    topo_gridx = collect(range(0, xsize, length=xnum))
    topo_gridy = zeros(xnum)
    add_basin_or_ridge!(
        topo_gridx, topo_gridy,
        topo_geometry.xo_topo,
        topo_geometry.width_topo, topo_geometry.height_topo
    )
    add_basin_or_ridge!(
        topo_gridx, topo_gridy,
        topo_geometry.xo_topo + topo_geometry.width_topo + topo_geometry.width_basin,
        topo_geometry.width_topo, topo_geometry.height_topo
    )
    gridx_b = collect(range(0, xsize, length=xnum))
    sediment_thickness_original = zeros(xnum)
    add_basin_or_ridge!(
        topo_gridx, sediment_thickness_original,
        topo_geometry.xo_topo + topo_geometry.width_topo,
        topo_geometry.width_basin, topo_geometry.sediment_thickness_original
    )
    return topo_gridx, topo_gridy, gridx_b, sediment_thickness_original
end

function add_basin_or_ridge!(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    feature_start::Float64,
    feature_width::Float64,
    feature_amp::Float64
)::Nothing
    dx_smooth = feature_width/5.0
    for i in 1:length(topo_gridx)
        x = topo_gridx[i]
        if feature_start < x < feature_start + feature_width
            if feature_start + dx_smooth <= x <= feature_start + feature_width - dx_smooth
                topo_gridy[i] += feature_amp
            elseif x < feature_start + dx_smooth
                topo_gridy[i] += feature_amp * (x - feature_start)/dx_smooth
            elseif x > feature_start + feature_width - dx_smooth
                topo_gridy[i] += feature_amp * 
                    (feature_start + feature_width - x)/dx_smooth
            end
        end
    end
    return nothing
end

function make_plot(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    top_sediment::Vector{Float64},
    bsmt_gridy::Vector{Float64};
    figsize::Tuple{Int,Int}=(10, 5)
)::Nothing
    dpi = 150
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)
    p = Plots.plot(
        topo_gridx, topo_gridy,
        label="Topography",
        color=:blue,
        size=figsize_pixels,
        margin=10Plots.mm,
        aspect_ratio=:auto,
        dpi=dpi,
        linestyle=:solid,
        linewidth=2.0,
        legend=:bottomright
    )
    Plots.plot!(
        p, topo_gridx, top_sediment, label="Top Sediment", color=:green,
        linestyle=:solid, linewidth=2.0
    )
    Plots.plot!(
        p, topo_gridx, bsmt_gridy, label="Initial Basement", color=:red,
        linestyle=:solid, linewidth=2.0
    )
    Plots.xlabel!("x (m)")
    Plots.ylabel!("Topography")
    Plots.title!("Multiple Eruptions")
    Plots.yaxis!(:flip)
    Plots.savefig(p, "multiple_eruptions.png")
    return nothing
end

function make_time_stamp(i::Int, timestep::Float64)::String
    model_time_yr = i * timestep
    model_time_my = model_time_yr/1e6
    return string(model_time_my) * " Myr"
end

end # module 