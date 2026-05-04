using EarthBox
using Test

@info "Julia thread/CPU info for this test run" Threads.nthreads() Sys.CPU_THREADS

# Test groups selected via the EARTHBOX_TEST_GROUP env var:
#   "fast" — unit tests + fast benchmarks (default; CI-friendly, < a few minutes)
#   "slow" — only the long benchmarks (Plasticity Kaus10, Elastic Slab; ~10+ min)
#   "all"  — everything
const TEST_GROUP = lowercase(get(ENV, "EARTHBOX_TEST_GROUP", "fast"))
const RUN_FAST = TEST_GROUP in ("fast", "all")
const RUN_SLOW = TEST_GROUP in ("slow", "all")

if !(TEST_GROUP in ("fast", "slow", "all"))
    error("EARTHBOX_TEST_GROUP must be one of \"fast\", \"slow\", \"all\"; got \"$TEST_GROUP\"")
end

@testset "EarthBox.jl" verbose=true begin
    if RUN_FAST
        include("tests/aqua_tests.jl")
        include("tests/path_validation_tests.jl")
        include("tests/runtools_tests.jl")
        include("tests/parameter_macros_tests.jl")
        include("tests/parameter_registry_guard_test.jl")
        include("tests/array_macros_tests.jl")
        include("tests/array_registry_guard_test.jl")
        include("tests/bisection_interp_tests.jl")
        include("tests/compaction_tests.jl")
        include("tests/marker_compaction_tests.jl")
        include("tests/rock_props_tests.jl")
        include("tests/half_space_cooling_tests.jl")
        include("tests/serpentinization_tests.jl")
        include("tests/melt_drainage_divides_tests.jl")
        include("tests/melt_model_tests.jl")
        include("tests/sediment_transport_tests.jl")
        include("tests/lava_flow_tests.jl")
        include("tests/base_level_tests.jl")
        include("tests/gravity_tests.jl")
        include("tests/benchmark_tests_fast.jl")
    end
    if RUN_SLOW
        include("tests/benchmark_tests_slow.jl")
    end
end
