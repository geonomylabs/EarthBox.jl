# Fast regression for the Stokes sinker multigrid path (3D, small grid).
#
# Benchmark / profiling (full-size StokesSinker):
#   ENV["EARTHBOX_MG_TIMING"] = "1"   # per-phase V-cycle time split (see MultigridVCycle)
#   ENV["EARTHBOX_MG_DEBUG"] = "0"    # omit final solution statistics printout
#
# Run this file:  julia --project=. test/multigrid/stokes_sinker_smoke.jl

using EarthBox
using Test

ENV["EARTHBOX_MG_DEBUG"] = "0"

import EarthBox.StokesContinuitySolver.MultigridManager.MultigridTests.StokesSinker as StokesSinker

@testset "Stokes sinker multigrid smoke (3D)" begin
    @test nothing === StokesSinker.run_stokes_sinker(smoke_test=true, model_type=:ThreeDimensional)
end
