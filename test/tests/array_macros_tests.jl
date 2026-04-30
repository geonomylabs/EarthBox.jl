using Test
using EarthBox

import EarthBox.Arrays.ArrayRegistry: ArrayData
import EarthBox.Arrays.ArrayMacros: @arrays
import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState

@testset "ArrayMacros.@arrays" begin

    @testset "Single ArrayData entry" begin
        nt = @arrays (
            gridx_b = ArrayData("m", GridArray1DState, "one-dimensional",
                                "x-coordinates of basic grid nodes"),
        )
        @test keys(nt) == (:gridx_b,)
        @test nt.gridx_b isa ArrayData
        @test nt.gridx_b.name        == "gridx_b"
        @test nt.gridx_b.units       == "m"
        @test nt.gridx_b.type        === GridArray1DState
        @test nt.gridx_b.grid_type   == "one-dimensional"
        @test nt.gridx_b.description == "x-coordinates of basic grid nodes"
    end

    @testset "Multiple entries" begin
        nt = @arrays (
            gridx_b = ArrayData("m",    GridArray1DState,   "one-dimensional", "x-coords"),
            ha0     = ArrayData("None", ScalarArray2DState, "basic",           "adiabatic backup"),
        )
        @test keys(nt) == (:gridx_b, :ha0)
        @test nt.gridx_b.name == "gridx_b"
        @test nt.ha0.name     == "ha0"
        @test nt.gridx_b.type === GridArray1DState
        @test nt.ha0.type     === ScalarArray2DState
        @test nt.gridx_b.grid_type == "one-dimensional"
        @test nt.ha0.grid_type     == "basic"
    end

    @testset "Multi-line star-concatenated description" begin
        nt = @arrays (
            big_array = ArrayData(
                "m", GridArray1DState, "one-dimensional",
                "First line of description "
                * "second line continues here "
                * "third line ends it.",
            ),
        )
        @test nt.big_array.name == "big_array"
        @test occursin("First line",     nt.big_array.description)
        @test occursin("second line",    nt.big_array.description)
        @test occursin("third line",     nt.big_array.description)
        @test nt.big_array.description ==
            "First line of description " *
            "second line continues here " *
            "third line ends it."
    end

    @testset "Field access still works" begin
        nt = @arrays (
            foo = ArrayData("m", GridArray1DState, "one-dimensional", "alpha"),
            bar = ArrayData("m", GridArray1DState, "one-dimensional", "beta"),
        )
        @test nt.foo.description == "alpha"
        @test nt.bar.description == "beta"
        for (k, v) in pairs(nt)
            @test String(k) == v.name
        end
    end

    @testset "Equivalence: @arrays vs hand-written NamedTuple" begin
        macro_nt = @arrays (
            gridx_b = ArrayData("m",    GridArray1DState,   "one-dimensional", "x-coords"),
            ha0     = ArrayData("None", ScalarArray2DState, "basic",           "adiabatic"),
        )
        hand_nt = (
            gridx_b = ArrayData("gridx_b", "m",    GridArray1DState,   "one-dimensional", "x-coords"),
            ha0     = ArrayData("ha0",     "None", ScalarArray2DState, "basic",           "adiabatic"),
        )
        @test keys(macro_nt) == keys(hand_nt)
        for k in keys(hand_nt)
            mp = macro_nt[k]
            hp = hand_nt[k]
            @test typeof(mp) === typeof(hp)
            @test mp.name        == hp.name
            @test mp.units       == hp.units
            @test mp.type        === hp.type
            @test mp.grid_type   == hp.grid_type
            @test mp.description == hp.description
        end
    end

    @testset "Negative: 5-arg call rejected" begin
        @test_throws LoadError @eval @arrays (
            a = ArrayData("a", "m", GridArray1DState, "one-dimensional", "desc"),
        )
    end

    @testset "Negative: unsupported constructor rejected" begin
        @test_throws LoadError @eval @arrays (
            a = SomeOtherThing("m", GridArray1DState, "one-dimensional", "desc"),
        )
    end

    @testset "Negative: non-tuple input rejected" begin
        @test_throws LoadError @eval @arrays 5
    end

    @testset "Negative: bare RHS without key rejected" begin
        @test_throws LoadError @eval @arrays (
            ArrayData("m", GridArray1DState, "one-dimensional", "desc"),
        )
    end

end
