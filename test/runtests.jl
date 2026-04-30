using EarthBox
using Test

@testset "EarthBox.jl" verbose=true begin
    include("tests/parameter_macros_tests.jl")
    include("tests/parameter_registry_guard_test.jl")
    include("tests/array_macros_tests.jl")
    include("tests/array_registry_guard_test.jl")
    include("tests/melt_drainage_divides_tests.jl")
    include("tests/benchmark_tests.jl")
end