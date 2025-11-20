using EarthBox
using Test

@testset "EarthBox.jl" begin
    expected_divides = [0.0, 18000.0, 84750.0, 150250.0, 215500.0, 281750.0, 348750.0, 416000.0, 483500.0, 500000.0]
    output_divides = TestManager.MeltDrainageDividesTest.run_test()

    @test length(output_divides) == length(expected_divides)
    for (out, expect) in zip(output_divides, expected_divides)
        @test isapprox(out, expect; atol=0.01)
    end

    (
        status_str, 
        max_relative_error_percentage, 
        relative_error_limit_percentage
    ) = BenchmarksManager.run_benchmark(:rayleigh_taylor_instability)
    @test status_str == "Success"
    @test max_relative_error_percentage < relative_error_limit_percentage


    (
        status_str, 
        max_relative_error_percentage, 
        relative_error_limit_percentage
    ) = BenchmarksManager.run_benchmark(:channel_flow_non_newtonian)
    @test status_str == "Success"
    @test max_relative_error_percentage < relative_error_limit_percentage

    (
        status_str, 
        max_relative_error_percentage, 
        relative_error_limit_percentage
    ) = BenchmarksManager.run_benchmark(:couette_flow_viscous_heating)
    @test status_str == "Success"
    @test max_relative_error_percentage < relative_error_limit_percentage

    (
        status_str, 
        max_relative_error_percentage, 
        relative_error_limit_percentage
    ) = BenchmarksManager.run_benchmark(:channel_flow_variable_conductivity)
    @test status_str == "Success"
    @test max_relative_error_percentage < relative_error_limit_percentage

    (
        status_str, 
        max_relative_error_percentage, 
        relative_error_limit_percentage
    ) = BenchmarksManager.run_benchmark(:solid_body_rotation)
    @test status_str == "Success"
    @test max_relative_error_percentage < relative_error_limit_percentage

end