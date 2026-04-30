using EarthBox
using Test

@testset "MeltDrainageDividesTest" begin
    expected_divides = [0.0, 18000.0, 84750.0, 150250.0, 215500.0, 281750.0, 
                        348750.0, 416000.0, 483500.0, 500000.0]
    output_divides = TestManager.MeltDrainageDividesTest.run_test()

    @test length(output_divides) == length(expected_divides)
    for (out, expect) in zip(output_divides, expected_divides)
        @test isapprox(out, expect; atol=0.01)
    end
end