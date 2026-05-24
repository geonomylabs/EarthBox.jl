using EarthBox
using Test

import EarthBox.MeltModel.Drainage: redistribute_extrusion_volumes!

""" Run the pure redistribution and return the new extrusion volumes.

`vols_o`, `xstart_o`, `xend_o` describe the old basins; `xstart`, `xend`
describe the new basins. Arrays are sized to their basin counts.
"""
function run_redistribution(vols_o, xstart_o, xend_o, xstart, xend)
    nnew = length(xstart)
    nold = length(xstart_o)
    vols = zeros(Float64, nnew)
    redistribute_extrusion_volumes!(
        vols, collect(Float64, vols_o),
        collect(Float64, xstart), collect(Float64, xend),
        collect(Float64, xstart_o), collect(Float64, xend_o),
        nnew, nold,
    )
    return vols
end

@testset "ExtrusionVolumeRedistributionTest" begin

    # Regression for the doubling bug: a full-width basin [0, 500000] (midpoint
    # 250000) splits at its own midpoint. The old midpoint lands exactly on the
    # shared boundary of the two new basins. It must be counted once (not twice)
    # and, by the right-favouring tie-break, go to the right basin.
    @testset "full-width split at midpoint is not doubled" begin
        V = 1.0e7
        vols = run_redistribution(
            [V], [0.0], [500_000.0],
            [0.0, 250_000.0], [250_000.0, 500_000.0],
        )
        @test isapprox(sum(vols), V)            # not 2V
        @test isapprox(vols[1], 0.0)            # left basin
        @test isapprox(vols[2], V)              # right basin (tie-break)
    end

    # Merge two basins back into one full-width basin: both volumes land in the
    # single new basin, total conserved.
    @testset "merge 2 -> 1 conserves total" begin
        V1, V2 = 7.0e6, 3.0e6
        vols = run_redistribution(
            [V1, V2], [0.0, 250_000.0], [250_000.0, 500_000.0],
            [0.0], [500_000.0],
        )
        @test isapprox(sum(vols), V1 + V2)
        @test isapprox(vols[1], V1 + V2)
    end

    # The 3 -> 2 transition with a thin sliver basin, taken from the case37 log
    # ([0,249800],[249800,250000],[250000,500000] -> [0,249800],[249800,500000]).
    @testset "3 -> 2 with thin sliver conserves total" begin
        a, b, c = 4.0e6, 1.0e5, 5.0e6
        vols = run_redistribution(
            [a, b, c],
            [0.0, 249_800.0, 250_000.0], [249_800.0, 250_000.0, 500_000.0],
            [0.0, 249_800.0], [249_800.0, 500_000.0],
        )
        @test isapprox(sum(vols), a + b + c)
        @test isapprox(vols[1], a)              # only midpoint 124900 -> left
        @test isapprox(vols[2], b + c)          # midpoints 249900, 375000 -> right
    end

    # Identity: same basins in and out leaves each basin's volume unchanged.
    @testset "identity transition preserves per-basin volume" begin
        v = [2.0e6, 4.0e6, 6.0e6]
        xs = [0.0, 150_000.0, 350_000.0]
        xe = [150_000.0, 350_000.0, 500_000.0]
        vols = run_redistribution(v, xs, xe, xs, xe)
        @test isapprox(sum(vols), sum(v))
        for k in 1:3
            @test isapprox(vols[k], v[k])
        end
    end

    # General non-coincident split (midpoints fall strictly inside new basins).
    @testset "non-coincident split conserves total" begin
        V = 8.0e6
        vols = run_redistribution(
            [V], [0.0], [500_000.0],
            [0.0, 300_000.0], [300_000.0, 500_000.0],
        )
        @test isapprox(sum(vols), V)
        @test isapprox(vols[1], V)              # midpoint 250000 in [0,300000)
        @test isapprox(vols[2], 0.0)
    end
end
