using Test
using EarthBox
using EarthBox.InputTools.InputBounds: validate_value, validate_total_markers,
                                       TOTAL_MARKERS_CAP

@testset "InputBounds" begin
    @testset "integer bounds — accepts in-range" begin
        validate_value("xnum", 301)
        validate_value("ynum", 61)
        validate_value("ntimestep_max", 550)
        validate_value("nglobal", 1)
        validate_value("nsand_layers", 24)
        validate_value("ntime_increase_1", 10)
        @test true  # reaching here means no throws
    end

    @testset "integer bounds — rejects out-of-range" begin
        @test_throws ArgumentError validate_value("xnum", 0)
        @test_throws ArgumentError validate_value("xnum", -1)
        @test_throws ArgumentError validate_value("xnum", 2_000_000)
        @test_throws ArgumentError validate_value("ynum", 10_000_000)
        @test_throws ArgumentError validate_value("ntimestep_max", -1)
        @test_throws ArgumentError validate_value("ntimestep_max", 100_000_000)
        @test_throws ArgumentError validate_value("nglobal", 100_000)
        @test_throws ArgumentError validate_value("nsand_layers", 5_000)
    end

    @testset "integer bounds — rejects non-integer" begin
        @test_throws ArgumentError validate_value("xnum", 301.5)
        @test_throws ArgumentError validate_value("ynum", "61")
    end

    @testset "real bounds — accepts in-range" begin
        validate_value("nmarkers_cell_x", 4.0)
        validate_value("nmarkers_cell_y", 4)
        validate_value("nmarkers_cell_x", 100)
        @test true
    end

    @testset "real bounds — rejects out-of-range" begin
        @test_throws ArgumentError validate_value("nmarkers_cell_x", 0.5)
        @test_throws ArgumentError validate_value("nmarkers_cell_x", 101.0)
        @test_throws ArgumentError validate_value("nmarkers_cell_y", -1.0)
    end

    @testset "real bounds — rejects NaN/Inf" begin
        @test_throws ArgumentError validate_value("nmarkers_cell_x", NaN)
        @test_throws ArgumentError validate_value("nmarkers_cell_x", Inf)
        @test_throws ArgumentError validate_value("nmarkers_cell_y", -Inf)
    end

    @testset "FLOAT_FINITE rejects NaN/Inf" begin
        @test_throws ArgumentError validate_value("xsize", NaN)
        @test_throws ArgumentError validate_value("ysize", Inf)
        @test_throws ArgumentError validate_value("viscosity_min", NaN)
        @test_throws ArgumentError validate_value("timestep_viscoelastic", -Inf)
        @test_throws ArgumentError validate_value("max_temp_change", NaN)
        @test_throws ArgumentError validate_value("gravity_y", Inf)
    end

    @testset "FLOAT_FINITE accepts finite reals (incl. zero, negative for non-positive fields)" begin
        validate_value("gravity_x", 0.0)
        validate_value("gravity_x", -1.0)
        validate_value("gravity_y", 9.8)
        validate_value("temperature_top", 273.0)
        validate_value("temperature_top", -50.0)  # not POSITIVE, just finite
        validate_value("pressure_bc", 20.0)
        @test true
    end

    @testset "FLOAT_POSITIVE rejects <= 0" begin
        @test_throws ArgumentError validate_value("xsize", 0.0)
        @test_throws ArgumentError validate_value("xsize", -0.1)
        @test_throws ArgumentError validate_value("ysize", 0.0)
        @test_throws ArgumentError validate_value("viscosity_min", 0.0)
        @test_throws ArgumentError validate_value("viscosity_min", -1.0)
        @test_throws ArgumentError validate_value("timestep_viscoelastic", 0.0)
        @test_throws ArgumentError validate_value("max_temp_change", -1.0)
    end

    @testset "FLOAT_POSITIVE accepts > 0" begin
        validate_value("xsize", 0.30)
        validate_value("ysize", 0.06)
        validate_value("viscosity_min", 1.0e4)
        validate_value("timestep_viscoelastic", 1.9e-9)
        validate_value("max_temp_change", 70.0)
        @test true
    end

    @testset "unknown names pass silently" begin
        validate_value("unknown_field_xyz", 12345)
        validate_value("another_unknown", "string value")
        validate_value("yet_another", NaN)
        @test true
    end

    @testset "validate_total_markers" begin
        # Sample (sandbox_extension example) — well under the cap.
        validate_total_markers(301, 61, 4.0, 4.0)

        # At cap boundary: 10^10 exactly should pass; cap+1 should fail.
        validate_total_markers(100_000, 100_000, 1, 1)  # 10^10 — passes
        @test_throws ArgumentError validate_total_markers(100_000, 100_000, 1, 2)
        @test_throws ArgumentError validate_total_markers(1_000_000, 1_000_000, 100, 100)
    end

    @testset "real example model.yml files pass validation" begin
        for example in [
            "examples/models/sandbox_extension/model.yml",
            "examples/models/slab_retreat/model.yml",
            "examples/models/sandbox_shortening/model.yml",
        ]
            result = EarthBox.InputTools.Reader.get_parameters_input_dict(example)
            @test result !== nothing
            @test length(result) > 0
        end
    end

    @testset "DoS-class config rejected at parse time" begin
        mktempdir(; prefix = "ib_dos_") do dir
            path = joinpath(dir, "evil.yml")
            open(path, "w") do io
                write(io, """
                    EarthBoxDimensions:
                      xnum: [5000000000, None, ' Way too large ']
                      ynum: [61, None, ' y ']
                      xsize: [0.30, m, ' width ']
                      ysize: [0.06, m, ' height ']
                      nmarkers_cell_x: [4.0, None, ' nmcx ']
                      nmarkers_cell_y: [4.0, None, ' nmcy ']
                    """)
            end
            @test_throws ArgumentError EarthBox.InputTools.Reader.get_parameters_input_dict(path)
        end
    end

    @testset "NaN in finite-required field rejected at parse time" begin
        mktempdir(; prefix = "ib_nan_") do dir
            path = joinpath(dir, "nan.yml")
            open(path, "w") do io
                write(io, """
                    EarthBoxDimensions:
                      xsize: [.nan, m, ' nan ']
                    """)
            end
            @test_throws ArgumentError EarthBox.InputTools.Reader.get_parameters_input_dict(path)
        end
    end

    @testset "total markers cap rejected via real reader path" begin
        mktempdir(; prefix = "ib_total_") do dir
            path = joinpath(dir, "total.yml")
            open(path, "w") do io
                write(io, """
                    EarthBoxDimensions:
                      xnum: [500000, None, ' x ']
                      ynum: [500000, None, ' y ']
                      nmarkers_cell_x: [10.0, None, ' nmcx ']
                      nmarkers_cell_y: [10.0, None, ' nmcy ']
                    """)
            end
            # 500_000 * 500_000 * 10 * 10 = 2.5e13, well over 10^10 cap.
            @test_throws ArgumentError EarthBox.InputTools.Reader.get_parameters_input_dict(path)
        end
    end
end
