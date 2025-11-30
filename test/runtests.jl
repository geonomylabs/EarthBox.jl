using EarthBox
using Test

# Check if big tests should be run
const RUN_BIG_TESTS = get(ENV, "EARTHBOX_RUN_BIG_TESTS", "false") == "true"

@testset "EarthBox.jl" begin
    @testset "Basic Tests" begin
        expected_divides = [0.0, 18000.0, 84750.0, 150250.0, 215500.0, 281750.0, 
                           348750.0, 416000.0, 483500.0, 500000.0]
        output_divides = TestManager.MeltDrainageDividesTest.run_test()

        @test length(output_divides) == length(expected_divides)
        for (out, expect) in zip(output_divides, expected_divides)
            @test isapprox(out, expect; atol=0.01)
        end
    end

    if RUN_BIG_TESTS
        @testset "Benchmarks" begin
            @testset "Rayleigh-Taylor Instability" begin
                (
                    status_str, 
                    max_relative_error_percentage, 
                    relative_error_limit_percentage
                ) = BenchmarksManager.run_benchmark(:rayleigh_taylor_instability)
                @test status_str == "Success"
                @test max_relative_error_percentage < relative_error_limit_percentage
            end

            @testset "Channel Flow Non-Newtonian" begin
                (
                    status_str, 
                    max_relative_error_percentage, 
                    relative_error_limit_percentage
                ) = BenchmarksManager.run_benchmark(:channel_flow_non_newtonian)
                @test status_str == "Success"
                @test max_relative_error_percentage < relative_error_limit_percentage
            end

            @testset "Couette Flow Viscous Heating" begin
                (
                    status_str, 
                    max_relative_error_percentage, 
                    relative_error_limit_percentage
                ) = BenchmarksManager.run_benchmark(:couette_flow_viscous_heating)
                @test status_str == "Success"
                @test max_relative_error_percentage < relative_error_limit_percentage
            end

            @testset "Channel Flow Variable Conductivity" begin
                (
                    status_str, 
                    max_relative_error_percentage, 
                    relative_error_limit_percentage
                ) = BenchmarksManager.run_benchmark(:channel_flow_variable_conductivity)
                @test status_str == "Success"
                @test max_relative_error_percentage < relative_error_limit_percentage
            end

            @testset "Solid Body Rotation" begin
                (
                    status_str, 
                    max_relative_error_percentage, 
                    relative_error_limit_percentage
                ) = BenchmarksManager.run_benchmark(:solid_body_rotation)
                @test status_str == "Success"
                @test max_relative_error_percentage < relative_error_limit_percentage
            end
        end
    else
        @info "Skipping big tests. Set EARTHBOX_RUN_BIG_TESTS=true to run them."
    end
end