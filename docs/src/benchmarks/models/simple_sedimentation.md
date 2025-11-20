# Simple Sedimentation

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:simple_sedimentation);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 1.92316
Max relative error: 0.060 %
Max relative error limit: 10.000 %

Test Successful: max relative error is not greater than limit: 0.060% <= 10.000%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_17-20-26/simple_sedimentation_output
```

![Simple Sedimentation Benchmark (Topography)](images/simple_sedimentation_topo_5.png)
##### fig:simple-sedimentation
*Benchmark result using EarthBox version 0.1.0*
