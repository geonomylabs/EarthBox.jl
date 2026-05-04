using Test
using EarthBox.PathValidation

@testset "PathValidation" begin
    @testset "accepts paths under default prefixes" begin
        @test PathValidation.validate_path("/home/user/foo") == "/home/user/foo"
        @test PathValidation.validate_path("/mnt/data/bar") == "/mnt/data/bar"
    end

    @testset "accepts filenames containing '..' substring" begin
        # Regression: substring '..' check was over-strict; legitimate filenames
        # with stacked dots must pass.
        @test PathValidation.validate_path("/home/user/report..final.txt") ==
              "/home/user/report..final.txt"
    end

    @testset "rejects directory traversal" begin
        @test_throws ArgumentError PathValidation.validate_path("../etc/passwd")
        @test_throws ArgumentError PathValidation.validate_path("/home/../etc/passwd")
        @test_throws ArgumentError PathValidation.validate_path("a/../../b")
    end

    @testset "rejects absolute paths outside allowed prefixes" begin
        @test_throws ArgumentError PathValidation.validate_path("/etc/passwd")
        @test_throws ArgumentError PathValidation.validate_path("/tmp/foo")
        @test_throws ArgumentError PathValidation.validate_path("/var/log/syslog")
    end

    @testset "configurable allowed_prefixes" begin
        @test PathValidation.validate_path("/tmp/foo";
                                           allowed_prefixes = ["/tmp"]) == "/tmp/foo"
        @test PathValidation.validate_path("/opt/julia";
                                           allowed_prefixes = ["/opt"]) == "/opt/julia"
        @test_throws ArgumentError PathValidation.validate_path(
            "/home/user/foo"; allowed_prefixes = ["/tmp"])
    end

    @testset "relative paths are not prefix-checked" begin
        @test PathValidation.validate_path("foo/bar") == "foo/bar"
        @test PathValidation.validate_path("./foo") == "foo"
    end

    @testset "resolve_symlinks rejects escapes via symlink" begin
        mktempdir(prefix = "pv_outside_") do outside_dir
            mktempdir(; prefix = "pv_inside_") do inside_dir
                # The inside_dir is under the system tmpdir; configure prefixes
                # so inside_dir is allowed but outside_dir is not.
                inside_prefix = inside_dir
                link_path = joinpath(inside_dir, "escape")
                symlink(outside_dir, link_path)

                # Without symlink resolution: passes (link path is under prefix).
                @test PathValidation.validate_path(link_path;
                                                   allowed_prefixes = [inside_prefix]) ==
                      link_path

                # With symlink resolution: rejected because realpath escapes.
                @test_throws ArgumentError PathValidation.validate_path(
                    link_path;
                    allowed_prefixes = [inside_prefix],
                    resolve_symlinks = true)
            end
        end
    end

    @testset "validate_directory_path requires existence" begin
        mktempdir(; prefix = "pv_dir_") do dir
            @test PathValidation.validate_directory_path(
                dir; allowed_prefixes = [dirname(dir)]) == dir
            missing_dir = joinpath(dir, "nope")
            @test_throws ArgumentError PathValidation.validate_directory_path(
                missing_dir; allowed_prefixes = [dirname(dir)])
        end
    end

    @testset "validate_file_path requires existence" begin
        mktempdir(; prefix = "pv_file_") do dir
            file = joinpath(dir, "ok.txt")
            write(file, "hello")
            @test PathValidation.validate_file_path(
                file; allowed_prefixes = [dirname(dir)]) == file
            @test_throws ArgumentError PathValidation.validate_file_path(
                joinpath(dir, "missing.txt"); allowed_prefixes = [dirname(dir)])
        end
    end
end
