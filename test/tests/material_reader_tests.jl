using Test
using EarthBox
import EarthBox.Markers.MarkerMaterials.MaterialsContainer.MaterialReader:
    validate_material_entry, read_materials_input

@testset "MaterialReader structural validation" begin
    # Helper: build a complete, valid raw entry (triplet-shaped values).
    function valid_entry()
        return Dict{Any, Any}(
            "mat_name"       => Any["SandboxStickyAir", "None", "name"],
            "mat_type"       => Any["General", "None", "type"],
            "mat_domain"     => Any["Atmosphere", "None", "domain"],
            "red_fraction"   => Any[0.0, "None", "r"],
            "green_fraction" => Any[1.0, "None", "g"],
            "blue_fraction"  => Any[1.0, "None", "b"],
        )
    end

    @testset "accepts a valid entry" begin
        validate_material_entry(1, valid_entry())
        @test true
    end

    @testset "rejects missing required key" begin
        for required_key in ["mat_name", "mat_type", "mat_domain",
                             "red_fraction", "green_fraction", "blue_fraction"]
            entry = valid_entry()
            delete!(entry, required_key)
            @test_throws ArgumentError validate_material_entry(2, entry)
        end
    end

    @testset "rejects malformed triplet (bare value)" begin
        entry = valid_entry()
        entry["mat_name"] = "SandboxStickyAir"
        @test_throws ArgumentError validate_material_entry(3, entry)
    end

    @testset "rejects triplet of wrong length" begin
        entry = valid_entry()
        entry["red_fraction"] = Any[0.0, "None"]  # missing description
        @test_throws ArgumentError validate_material_entry(4, entry)
    end

    @testset "rejects out-of-[0,1] RGB" begin
        for rgb_key in ["red_fraction", "green_fraction", "blue_fraction"]
            entry = valid_entry()
            entry[rgb_key] = Any[1.5, "None", "out of range"]
            @test_throws ArgumentError validate_material_entry(5, entry)

            entry = valid_entry()
            entry[rgb_key] = Any[-0.1, "None", "negative"]
            @test_throws ArgumentError validate_material_entry(5, entry)

            entry = valid_entry()
            entry[rgb_key] = Any[NaN, "None", "nan"]
            @test_throws ArgumentError validate_material_entry(5, entry)
        end
    end

    @testset "rejects non-real RGB" begin
        entry = valid_entry()
        entry["red_fraction"] = Any["0.5", "None", "string instead of float"]
        @test_throws ArgumentError validate_material_entry(6, entry)
    end

    @testset "rejects empty mat_name" begin
        entry = valid_entry()
        entry["mat_name"] = Any["", "None", "empty"]
        @test_throws ArgumentError validate_material_entry(7, entry)
    end

    @testset "rejects non-string mat_name" begin
        entry = valid_entry()
        entry["mat_name"] = Any[42, "None", "non-string"]
        @test_throws ArgumentError validate_material_entry(8, entry)
    end

    @testset "real example materials.yml passes" begin
        for example in [
            "examples/models/sandbox_extension/materials.yml",
            "examples/models/slab_retreat/materials.yml",
            "examples/models/sandbox_shortening/materials.yml",
        ]
            materials_dict, nmats = read_materials_input(example)
            @test nmats > 0
            # mat_name should be propagated as the first parameter on every entry.
            for (matid, params) in materials_dict
                @test haskey(params, "mat_name")
            end
        end
    end
end
