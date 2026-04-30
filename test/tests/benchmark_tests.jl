using EarthBox
using Test

@testset "BenchmarkTests" verbose=true begin
    delete_output = true

    @testset "Rayleigh-Taylor Instability" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :rayleigh_taylor_instability; delete_output=delete_output
        )
        @test status_str == "Success"
        @test max_relative_error_percentage < relative_error_limit_percentage
    end

    @testset "Channel Flow Non-Newtonian" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :channel_flow_non_newtonian; delete_output=delete_output
        )
        @test status_str == "Success"
        @test max_relative_error_percentage < relative_error_limit_percentage
    end

    @testset "Couette Flow Viscous Heating" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :couette_flow_viscous_heating; delete_output=delete_output
        )
        @test status_str == "Success"
        @test max_relative_error_percentage < relative_error_limit_percentage
    end

    @testset "Channel Flow Variable Conductivity" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :channel_flow_variable_conductivity; delete_output=delete_output
        )
        @test status_str == "Success"
        @test max_relative_error_percentage < relative_error_limit_percentage
    end

    @testset "Solid Body Rotation" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :solid_body_rotation; delete_output=delete_output
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

    @testset "Simple Sedimentation" begin
        (
            status_str,
            max_relative_error_percentage,
            relative_error_limit_percentage
        ) = BenchmarksManager.run_benchmark(
            :simple_sedimentation; delete_output=delete_output
        )
        @test status_str == "Success"
        @test max_relative_error_percentage < relative_error_limit_percentage
    end


end