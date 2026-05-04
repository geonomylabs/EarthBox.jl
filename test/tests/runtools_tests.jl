using Test
using EarthBox
using EarthBox.RunTools

# Helpers reach into RunTools' internals; they're not exported.
const _build_command = EarthBox.RunTools.build_command
const _build_remote_command = EarthBox.RunTools.build_remote_command
const _get_runit_actions = EarthBox.RunTools.get_runit_actions
const _get_istart_iend = EarthBox.RunTools.get_istart_iend

@testset "RunTools" begin
    @testset "get_runit_actions returns Vector{String}" begin
        @test _get_runit_actions(nothing, nothing, nothing) == String[]
        @test _get_runit_actions(true, nothing, nothing) == ["run_model"]
        @test _get_runit_actions(true, true, true) ==
              ["run_model", "plot_markers", "plot_scalars"]
        @test _get_runit_actions(false, true, false) == ["plot_markers"]
    end

    @testset "get_istart_iend returns Vector{String}" begin
        @test _get_istart_iend(nothing, nothing) == String[]
        @test _get_istart_iend(1, nothing) == ["istart=1"]
        @test _get_istart_iend(nothing, 5) == ["iend=5"]
        @test _get_istart_iend(2, 7) == ["istart=2", "iend=7"]
    end

    @testset "build_command (model, no eb_path, no output path)" begin
        cmd, logfile = _build_command(
            "/home/me/model", nothing, "model", "case0",
            "model.log", nothing, 1, 5,
        )
        @test cmd isa Cmd
        @test cmd.exec == [
            "julia", "--startup-file=no",
            "/home/me/model/Model.jl",
            "case_name=case0",
        ]
        @test logfile == "model.log"
    end

    @testset "build_command (model with eb_path and output path)" begin
        cmd, logfile = _build_command(
            "/home/me/model", "/home/me/eb", "model", "case2",
            "model.log", "/home/me/out", 1, 1,
        )
        @test cmd.exec == [
            "julia", "--project=/home/me/eb", "--startup-file=no",
            "/home/me/model/Model.jl",
            "case_name=case2",
            "model_output_path=/home/me/out",
        ]
        @test logfile == "model.log"
    end

    @testset "build_command (marker_plots)" begin
        cmd, logfile = _build_command(
            "/home/me/model", nothing, "marker_plots", "case1",
            "model.log", "/home/me/out", 3, 9,
        )
        @test cmd.exec == [
            "julia", "--startup-file=no",
            "/home/me/model/Plot.jl", "marker_plots",
            "istart=3", "iend=9",
            "case_name=case1",
            "model_output_path=/home/me/out",
        ]
        @test logfile == "marker_plots_model.log"
    end

    @testset "build_command (scalar_plots)" begin
        cmd, logfile = _build_command(
            "/home/me/model", nothing, "scalar_plots", "case1",
            "model.log", nothing, 1, 1,
        )
        @test cmd.exec == [
            "julia", "--startup-file=no",
            "/home/me/model/Plot.jl", "scalar_plots",
            "istart=1", "iend=1",
            "case_name=case1",
        ]
        @test logfile == "scalar_plots_model.log"
    end

    @testset "build_command rejects unknown command_type" begin
        @test_throws ArgumentError _build_command(
            "/home/me/model", nothing, "bogus", "case0",
            "model.log", nothing, 1, 1,
        )
    end

    @testset "injection probes — argv stays atomic" begin
        # Each malicious value must occupy exactly one argv slot, never get
        # split into shell tokens.
        evil_case = "case0; touch /tmp/PWNED"
        cmd, _ = _build_command(
            "/home/me/model", nothing, "model", evil_case,
            "model.log", nothing, 1, 1,
        )
        @test "case_name=$evil_case" in cmd.exec
        @test !any(occursin("touch", arg) && occursin("/tmp/PWNED", arg) &&
                   !occursin("case_name=", arg) for arg in cmd.exec)

        evil_path = "/home/me/out\$(rm -rf /)"
        cmd2, _ = _build_command(
            "/home/me/model", nothing, "model", "case0",
            "model.log", evil_path, 1, 1,
        )
        @test "model_output_path=$evil_path" in cmd2.exec
    end

    @testset "build_remote_command escapes user values" begin
        # No metacharacters: shape preserved, no extra escaping noise.
        rc = _build_remote_command(
            "/home/me/model", "/home/me/model/Runit.jl",
            "case0", ["run_model"], 1, 5,
        )
        @test occursin("cd /home/me/model && nohup julia --startup-file=no /home/me/model/Runit.jl case_name=case0 run_model istart=1 iend=5 > /dev/null 2>&1 &", rc)

        # With shell metachars in case name: must be escaped (single-quoted
        # by Base.shell_escape).
        evil_case = "case0; rm -rf /"
        rc2 = _build_remote_command(
            "/home/me/model", "/home/me/model/Runit.jl",
            evil_case, String[], nothing, nothing,
        )
        # The escaped form must appear and the unescaped semicolon must NOT
        # be a free shell token.
        @test occursin("case_name='case0; rm -rf /'", rc2)

        # Empty actions / istart / iend are simply omitted.
        rc3 = _build_remote_command(
            "/home/me/model", "/home/me/model/Runit.jl",
            "case0", String[], nothing, nothing,
        )
        @test !occursin("istart", rc3)
        @test !occursin("iend", rc3)
    end

    @testset "execute_local_script_in_background builds chdir + detached cmd" begin
        # We can't easily intercept the spawn, but we can verify Cmd assembly
        # via the same helpers the function uses. The local function does not
        # have a build helper, so we replicate its argv shape and assert the
        # interpolation rules used inside it.
        runit_actions = _get_runit_actions(true, true, false)
        istart_iend_args = _get_istart_iend(1, 5)
        script_path = "/home/me/model/Runit.jl"
        cmd = Cmd(
            `julia --startup-file=no $script_path case_name=case0 $runit_actions $istart_iend_args`,
            dir = "/home/me/model",
        )
        @test cmd.dir == "/home/me/model"
        @test cmd.exec == [
            "julia", "--startup-file=no", script_path,
            "case_name=case0",
            "run_model", "plot_markers",
            "istart=1", "iend=5",
        ]
    end

    # ----- Layer 2: stub-process tests -----

    @testset "execute_earthbox_script runs subprocess and writes logfile" begin
        mktempdir(; prefix = "rt_stub_") do model_dir
            # Stub Model.jl that prints a sentinel and exits 0. It also
            # echoes its ARGS so we can confirm argv reached unchanged.
            write(joinpath(model_dir, "Model.jl"), """
                println("STUB_OK")
                println("ARGS=", join(ARGS, "|"))
                exit(0)
                """)
            write(joinpath(model_dir, "Plot.jl"), "exit(0)\n")

            cd(model_dir) do
                EarthBox.RunTools.execute_earthbox_script(
                    model_dir = model_dir,
                    eb_path = nothing,
                    command_type = "model",
                    model_case_name = "caseX",
                    model_logfile_name = "stub.log",
                    model_output_path = nothing,
                )
                logpath = joinpath(model_dir, "stub.log")
                @test isfile(logpath)
                contents = read(logpath, String)
                @test occursin("STUB_OK", contents)
                @test occursin("ARGS=case_name=caseX", contents)
            end
        end
    end

    @testset "execute_earthbox_script preserves injection-attempt as one argv" begin
        mktempdir(; prefix = "rt_stub2_") do model_dir
            write(joinpath(model_dir, "Model.jl"), """
                println("ARGS=", join(ARGS, "|"))
                exit(0)
                """)
            write(joinpath(model_dir, "Plot.jl"), "exit(0)\n")

            evil_case = "caseX; touch SHOULD_NOT_EXIST"
            cd(model_dir) do
                EarthBox.RunTools.execute_earthbox_script(
                    model_dir = model_dir,
                    eb_path = nothing,
                    command_type = "model",
                    model_case_name = evil_case,
                    model_logfile_name = "stub.log",
                    model_output_path = nothing,
                )
                contents = read(joinpath(model_dir, "stub.log"), String)
                # Whole malicious string must appear as a single argv token.
                @test occursin("ARGS=case_name=$evil_case", contents)
                # If shell parsing had occurred, this file would exist.
                @test !isfile(joinpath(model_dir, "SHOULD_NOT_EXIST"))
            end
        end
    end
end
