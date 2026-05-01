using EarthBox
using Test

import EarthBox.Compaction.CompactionTools: compact_or_decompact

@testset "Compaction (Athy compaction kernel)" begin
    porosity_initial = 0.4              # fraction
    depth_decay_term = 1.0 / 2500.0     # 1/m

    # Decompact a 620 m layer at 10 km depth back to the surface.
    top_initial    = 10_000.0
    bottom_initial = 10_620.0
    top_new        =      0.0
    thickness_new = compact_or_decompact(
        porosity_initial, depth_decay_term,
        top_initial, bottom_initial, top_new,
    )
    expected_thickness_decompact = 924.7998193672524

    @test isapprox(thickness_new, expected_thickness_decompact; atol=1e-6)

    # Round-trip: re-compact the decompacted layer back down to depth.
    top_round_trip = top_initial
    thickness_round_trip = compact_or_decompact(
        porosity_initial, depth_decay_term,
        top_new, top_new + thickness_new, top_round_trip,
    )
    expected_thickness_round_trip = 619.6100794442791

    @test isapprox(thickness_round_trip, expected_thickness_round_trip; atol=1e-6)
end
