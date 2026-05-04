using Test
using EarthBox.PathValidation: validate_safe_output_path, SYSTEM_DENY_PREFIXES

@testset "validate_safe_output_path" begin
    @testset "accepts common output roots" begin
        @test validate_safe_output_path("/home/me/runs") == "/home/me/runs"
        @test validate_safe_output_path("/Users/me/runs") == "/Users/me/runs"
        @test validate_safe_output_path("/mnt/data/runs") == "/mnt/data/runs"
        @test validate_safe_output_path("/scratch/proj") == "/scratch/proj"
        @test validate_safe_output_path("/data/case0_output") == "/data/case0_output"
        @test validate_safe_output_path("/opt/storage/runs") == "/opt/storage/runs"
        @test validate_safe_output_path("/tmp/case0") == "/tmp/case0"
        @test validate_safe_output_path("/var/tmp/runs") == "/var/tmp/runs"
    end

    @testset "accepts relative paths (no prefix check)" begin
        @test validate_safe_output_path("runs") == "runs"
        @test validate_safe_output_path("./runs") == "runs"
        @test validate_safe_output_path("output/case0") == "output/case0"
        @test validate_safe_output_path("subdir/x") == "subdir/x"
    end

    @testset "normalises paths" begin
        # Julia's normpath collapses ./ and duplicate slashes but preserves a
        # trailing slash; both forms resolve to the same directory either way.
        @test validate_safe_output_path("/home/me/./runs") == "/home/me/runs"
        @test validate_safe_output_path("/home//me/runs") == "/home/me/runs"
        @test validate_safe_output_path("/home/me/runs/") == "/home/me/runs/"
    end

    @testset "rejects directory traversal" begin
        @test_throws ArgumentError validate_safe_output_path("../etc/passwd")
        @test_throws ArgumentError validate_safe_output_path("../../foo")
        @test_throws ArgumentError validate_safe_output_path("a/../../b")
    end

    @testset "rejects absolute paths under system directories" begin
        for prefix in SYSTEM_DENY_PREFIXES
            @test_throws ArgumentError validate_safe_output_path(prefix)
            @test_throws ArgumentError validate_safe_output_path(prefix * "/foo")
            @test_throws ArgumentError validate_safe_output_path(prefix * "/sub/dir")
        end
    end

    @testset "rejects post-normalisation escapes into system dirs" begin
        @test_throws ArgumentError validate_safe_output_path("/home/me/../../etc/passwd")
        @test_throws ArgumentError validate_safe_output_path("/scratch/../proc/1")
    end

    @testset "does not match prefixes by substring" begin
        # Paths that look like system prefixes but aren't (e.g. /etcetera, /ushare)
        # must be accepted — only the exact prefix segment counts.
        @test validate_safe_output_path("/etcetera/runs") == "/etcetera/runs"
        @test validate_safe_output_path("/usher/runs") == "/usher/runs"
        @test validate_safe_output_path("/lib2/runs") == "/lib2/runs"
        @test validate_safe_output_path("/devops/runs") == "/devops/runs"
    end

    @testset "rejects empty path" begin
        @test_throws ArgumentError validate_safe_output_path("")
    end

    @testset "accepts AbstractString (e.g. SubString)" begin
        # CLI-arg parsers return SubString{String} after strip(); confirm it
        # passes through cleanly.
        s = strip("  /home/me/runs  ")
        @test validate_safe_output_path(s) == "/home/me/runs"
        @test validate_safe_output_path(s) isa String
    end
end
