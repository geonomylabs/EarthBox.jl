using Test
using EarthBox

import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParamMacros: @params

@testset "ParamMacros.@params" begin

    @testset "Single ParameterInt entry" begin
        nt = @params (a = ParameterInt(7, "None", "an int parameter"),)
        @test keys(nt) == (:a,)
        @test nt.a isa ParameterInt
        @test nt.a.name == "a"
        @test nt.a.value == 7
        @test nt.a.units == "None"
        @test nt.a.description == "an int parameter"
    end

    @testset "Single ParameterFloat entry" begin
        nt = @params (b = ParameterFloat(2.5, "m", "a float parameter"),)
        @test keys(nt) == (:b,)
        @test nt.b isa ParameterFloat
        @test nt.b.name == "b"
        @test nt.b.value == 2.5
        @test nt.b.units == "m"
        @test nt.b.description == "a float parameter"
    end

    @testset "Single ParameterStr entry" begin
        nt = @params (c = ParameterStr("hello", "None", "a string parameter"),)
        @test keys(nt) == (:c,)
        @test nt.c isa ParameterStr
        @test nt.c.name == "c"
        @test nt.c.value == "hello"
        @test nt.c.units == "None"
        @test nt.c.description == "a string parameter"
    end

    @testset "Mixed tuple with all three types" begin
        nt = @params (
            iuse_thing      = ParameterInt(0,    "None",   "flag: 0=off, 1=on"),
            thickness_m     = ParameterFloat(1.0e3, "m",   "layer thickness"),
            material_label  = ParameterStr("granite", "None", "material name"),
        )
        @test keys(nt) == (:iuse_thing, :thickness_m, :material_label)
        @test nt.iuse_thing.name     == "iuse_thing"
        @test nt.thickness_m.name    == "thickness_m"
        @test nt.material_label.name == "material_label"
        @test nt.iuse_thing isa ParameterInt
        @test nt.thickness_m isa ParameterFloat
        @test nt.material_label isa ParameterStr
        @test nt.iuse_thing.value     == 0
        @test nt.thickness_m.value    == 1000.0
        @test nt.material_label.value == "granite"
    end

    @testset "Multi-line star-concatenated description" begin
        nt = @params (
            plate_thickness = ParameterFloat(
                100000.0, "m",
                "Thickness of plate in meters used for inflow-outflow depth-dependent extension and "
                * "contraction boundary conditions"
            ),
        )
        @test nt.plate_thickness.name == "plate_thickness"
        @test occursin("inflow-outflow", nt.plate_thickness.description)
        @test occursin("contraction boundary conditions", nt.plate_thickness.description)
        @test nt.plate_thickness.description ==
            "Thickness of plate in meters used for inflow-outflow depth-dependent extension and " *
            "contraction boundary conditions"
    end

    @testset "Expression default value" begin
        nt = @params (
            sec_per_yr = ParameterFloat(365.25*24.0*3600.0, "s/yr", "Seconds per year"),
        )
        @test nt.sec_per_yr.name == "sec_per_yr"
        @test nt.sec_per_yr.value == 365.25 * 24.0 * 3600.0
    end

    @testset "Equivalence: @params vs hand-written NamedTuple" begin
        macro_nt = @params (
            iuse_thing  = ParameterInt(0,   "None", "flag"),
            thickness_m = ParameterFloat(1.0e3, "m", "thickness"),
            label       = ParameterStr("g", "None", "label"),
        )
        hand_nt = (
            iuse_thing  = ParameterInt(0, "iuse_thing", "None", "flag"),
            thickness_m = ParameterFloat(1.0e3, "thickness_m", "m", "thickness"),
            label       = ParameterStr("g", "label", "None", "label"),
        )
        @test keys(macro_nt) == keys(hand_nt)
        for k in keys(hand_nt)
            mp = macro_nt[k]
            hp = hand_nt[k]
            @test typeof(mp) === typeof(hp)
            @test mp.name        == hp.name
            @test mp.value       == hp.value
            @test mp.units       == hp.units
            @test mp.description == hp.description
        end
    end

    @testset "Field access still works" begin
        nt = @params (foo = ParameterInt(42, "None", "answer"),)
        @test nt.foo.value == 42
        @test nt.foo.name  == "foo"
        for (k, v) in pairs(nt)
            @test String(k) == v.name
        end
    end

    @testset "Negative: 4-arg call rejected" begin
        @test_throws LoadError @eval @params (
            a = ParameterInt(0, "a", "None", "desc"),
        )
    end

    @testset "Negative: unsupported constructor rejected" begin
        @test_throws LoadError @eval @params (
            a = SomeOtherThing(0, "None", "desc"),
        )
    end

    @testset "Negative: non-tuple input rejected" begin
        @test_throws LoadError @eval @params 5
    end

    @testset "Negative: bare RHS without key rejected" begin
        @test_throws LoadError @eval @params (
            ParameterInt(0, "None", "desc"),
        )
    end

end
