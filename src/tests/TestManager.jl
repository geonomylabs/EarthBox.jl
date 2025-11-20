"""
TestManager module contains all the tests for the EarthBox model.

Exported Test Modules:
- SedimentTransportTest
- LavaFlowTestSingle
- LavaFlowTestMultiple
- MarkerCompactionTest
- MagmaBodyTest
- MeltDrainageDividesTest
- MeltDamageTest
- GabbroSolidusLiquidusTest
- PeridotiteSolidusLiquidusTest
- ConductivityTest
- HeatCapacityTest
- SerpentinizationTest
- BaseLevelTest
- LatentHeatTest
- HalfSpaceCoolingTest
- BisectionInterpTest
- CompactionTest
- TestGridGravity
- TestCellGravity

To run tests, use the following syntax:

TestManager.[test_module].run_test()

For example,

```julia
using EarthBox
TestManager.SedimentTransportTest.run_test()
```

"""
module TestManager

include("tests/surface_processes/sediment_transport/SedimentTransportTest.jl")
include("tests/surface_processes/lava_flow/LavaFlowTestSingle.jl")
include("tests/surface_processes/lava_flow/LavaFlowTestMultiple.jl")
include("tests/compaction/marker_compaction/MarkerCompactionTest.jl")
include("tests/melt_model/extraction/MagmaBodyTest.jl")
include("tests/melt_model/divides/MeltDrainageDividesTest.jl")
include("tests/melt_model/melt_damage/MeltDamageTest.jl")
include("tests/melt_model/solidus_liquidus/GabbroSolidusLiquidusTest.jl")
include("tests/melt_model/solidus_liquidus/PeridotiteSolidusLiquidusTest.jl")
include("tests/rock_props/conductivity/ConductivityTest.jl")
include("tests/rock_props/heat_capacity/HeatCapacityTest.jl")
include("tests/serpentinization/SerpentinizationTest.jl")
include("tests/surface_processes/sealevel/base_level/BaseLevelTest.jl")
include("tests/melt_model/latent_heat/LatentHeatTest.jl")
include("tests/marker_temperature/half_space_cooling/HalfSpaceCoolingTest.jl")
include("tests/utils/bisection_interp/BisectionInterpTest.jl")
include("tests/compaction/CompactionTest.jl")
include("tests/plot_tools/ScatterTest.jl")
include("tests/gravity/TestCellGravity.jl")
include("tests/gravity/TestGridGravity.jl")

import .SedimentTransportTest, .LavaFlowTestSingle, .LavaFlowTestMultiple
import .MarkerCompactionTest, .MagmaBodyTest, .MeltDrainageDividesTest
import .MeltDamageTest, .GabbroSolidusLiquidusTest, .PeridotiteSolidusLiquidusTest
import .ConductivityTest, .HeatCapacityTest, .SerpentinizationTest
import .BaseLevelTest, .LatentHeatTest, .HalfSpaceCoolingTest
import .BisectionInterpTest, .CompactionTest, .TestGridGravity, .TestCellGravity

export SedimentTransportTest, LavaFlowTestSingle, LavaFlowTestMultiple
export MarkerCompactionTest, MagmaBodyTest, MeltDrainageDividesTest
export MeltDamageTest, GabbroSolidusLiquidusTest, PeridotiteSolidusLiquidusTest
export ConductivityTest, HeatCapacityTest, SerpentinizationTest
export BaseLevelTest, LatentHeatTest, HalfSpaceCoolingTest
export BisectionInterpTest, CompactionTest, TestGridGravity, TestCellGravity

end # module