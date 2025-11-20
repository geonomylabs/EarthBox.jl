module CompactionTest

import EarthBox.Compaction.CompactionTools: compact_or_decompact

function run_test()::Nothing
    # 1436.9516319947998
    porosity_initial = 0.4  # fraction
    depth_decay_term = 1.0/2500.0  # 1/m
    top_initial = 10_000.0  # meters
    bottom_initial = 10_620.0
    top_new = 0.0
    thickness_new = compact_or_decompact(
        porosity_initial,
        depth_decay_term,
        top_initial,
        bottom_initial,
        top_new
    )
    println("Thickness of compacted layer: ", thickness_new)

    return nothing
end

end # module 