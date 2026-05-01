using EarthBox
using Test

import EarthBox.MathTools: linear_interp_bisection

@testset "BisectionInterp" begin
    gridx = collect(range(0.0, 500_000.0, length=100))
    gridy = collect(range(3.0, 10.0, length=length(gridx)))

    sample_xs = [
        0.0,           # at left edge
        100_000.0,     # interior
        250_000.0,     # interior (middle)
        500_000.0,     # at right edge
        800_000.0,     # extrapolation right
        -50_000.0,     # extrapolation left
    ]
    expected_ys = [
        3.0,
        4.4,
        6.5,
        10.0,
        14.200000000000017,    # extrapolation right
        2.300000000000002,     # extrapolation left
    ]
    for (k, x) in enumerate(sample_xs)
        y = linear_interp_bisection(gridx, gridy, x)
        @test isapprox(y, expected_ys[k]; atol=1e-9)
    end
end
