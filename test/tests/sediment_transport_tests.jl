using EarthBox
using Test

import EarthBox.ConversionFuncs: years_to_seconds
import EarthBox.ConversionFuncs: meters_per_year_to_meters_per_seconds as m_yr_to_m_s
import EarthBox.ConversionFuncs: meters_squared_per_year_to_meters_squared_per_second as m2_yr_to_m2_s
import EarthBox.ConversionFuncs: mm_per_yr_to_meters_per_seconds as mm_yr_to_m_s
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.SurfaceProcesses.SedimentTransport: SedimentTransportSolverManager

function build_sediment_transport_solver()
    xsize = 500_000.0
    dx = 500.0
    xnum = floor(Int, xsize/dx) + 1
    @assert xnum == 1001

    xo_topo = 200_000.0
    width_topo = 100_000.0
    height_topo = 1000.0
    xo_basin = xo_topo + width_topo
    width_basin = 100_000.0
    xf_basin = xo_basin + width_basin
    depth_basin = 2000.0
    sediment_thickness_original = 2000.0
    y_sealevel = 0.0

    topo_gridx = collect(range(0, xsize, length=xnum))
    topo_gridy = zeros(xnum)
    sediment_thickness = zeros(xnum)
    for i in eachindex(topo_gridx)
        x = topo_gridx[i]
        if xo_topo <= x <= xo_topo + width_topo
            topo_gridy[i] -= height_topo
        end
        if xo_basin <= x <= xf_basin
            topo_gridy[i] += depth_basin
            sediment_thickness[i] = sediment_thickness_original
        end
    end

    gridx_b = collect(range(0, xsize, length=xnum))

    sediment_transport_parameters = SedimentTransportParameters(
        m2_yr_to_m2_s(0.25),
        m_yr_to_m_s(1.0),
        1e-4,
        m2_yr_to_m2_s(1e2),
        2000.0,
        years_to_seconds(1_000.0),
        years_to_seconds(5.0e6),
        0.4,
        1/2000.0,
    )

    pelagic_sedimentation_rate = mm_yr_to_m_s(0.0)

    return SedimentTransportSolverManager.SedimentTransportSolver(
        (gridx_b[1], gridx_b[end]),
        topo_gridx,
        topo_gridy,
        sediment_thickness,
        sediment_transport_parameters,
        y_sealevel,
        pelagic_sedimentation_rate;
        use_collections=false,
        use_print_debug=false,
        use_constant_diffusivity=false,
        use_compaction_correction=true,
        use_optimized_solver=true,
    )
end

@testset "SedimentTransportTest" begin
    solver = build_sediment_transport_solver()
    SedimentTransportSolverManager.run_sediment_transport_time_steps!(solver)

    @test length(solver.topo_gridy) == 1001
    @test length(solver.sediment_thickness_initial_compacted) == 1001
    @test length(solver.compaction_displacement_max) == 1001

    sample_indices = [1, 501, 651, 701, 751, 1001]

    expected_topo_gridy = [
        0.0,
        -999.989456572056,
        1472.812650035865,
        1927.8608519598274,
        1651.3854094642174,
        2.0979457689966234,
    ]
    expected_sediment_thickness_initial_compacted = [
        0.0,
        0.0,
        1833.8213684155476,
        1970.9644443764814,
        1879.9644401836736,
        0.0,
    ]
    expected_compaction_displacement_max = [
        0.0,
        0.0,
        166.17863158445243,
        29.035555623518576,
        120.03555981632644,
        0.0,
    ]

    for (k, idx) in enumerate(sample_indices)
        @test isapprox(solver.topo_gridy[idx],
                       expected_topo_gridy[k]; atol=1e-6)
        @test isapprox(solver.sediment_thickness_initial_compacted[idx],
                       expected_sediment_thickness_initial_compacted[k]; atol=1e-6)
        @test isapprox(solver.compaction_displacement_max[idx],
                       expected_compaction_displacement_max[k]; atol=1e-6)
    end
end
