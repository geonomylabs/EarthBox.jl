using EarthBox
using Test

# Long-running benchmarks (multiple minutes each, up to ~10 min for Elastic Slab).
# Skipped unless EARTHBOX_TEST_GROUP is "slow" or "all".
@testset "BenchmarkTests (slow)" verbose=true begin
    delete_output = true

    @testset "Plasticity Benchmark Kaus10" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :plasticity_benchmark_kaus10; delete_output=delete_output
        )
        @test status_str == "Success"
        @test max_relative_error_percentage < relative_error_limit_percentage
    end

    @testset "Viscoelastic Stress Buildup" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :viscoelastic_stress_buildup; delete_output=delete_output
        )
        @test status_str == "Success"
        @test max_relative_error_percentage < relative_error_limit_percentage
    end

    @testset "Elastic Slab" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :elastic_slab; delete_output=delete_output
        )
        @test status_str == "Success"
        @test max_relative_error_percentage < relative_error_limit_percentage
    end
end
