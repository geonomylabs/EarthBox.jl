using EarthBox
using Test

import EarthBox.MeltModel.Drainage: calculate_drainage_divides

function build_drainage_divides_topography()
    xsize = 500_000.0
    dx = 500.0
    amplitude = 5_000.0
    wavelength = 100_000.0
    flat_threshold = 2_500.0

    xnum = floor(Int, xsize/dx) + 1
    topo_gridx = collect(range(0, xsize, length=xnum))
    sine_values = amplitude .* sin.(topo_gridx ./ wavelength) .+
                  amplitude .* sin.(3π .* topo_gridx ./ wavelength)
    sine_values[sine_values .>  flat_threshold]        .=  flat_threshold
    sine_values[sine_values .< -flat_threshold * 2.0]  .= -flat_threshold * 2.0
    return topo_gridx, sine_values
end

@testset "MeltDrainageDividesTest" begin
    expected_divides = [0.0, 18000.0, 84750.0, 150250.0, 215500.0, 281750.0,
                        348750.0, 416000.0, 483500.0, 500000.0]

    topo_gridx, topo_gridy = build_drainage_divides_topography()
    output_divides = calculate_drainage_divides(topo_gridx, topo_gridy)

    @test length(output_divides) == length(expected_divides)
    for (out, expect) in zip(output_divides, expected_divides)
        @test isapprox(out, expect; atol=0.01)
    end
end
