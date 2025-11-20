# Channel Flow Non-Newtonian

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:channel_flow_non_newtonian);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 0.08735
Max relative error: 1.008 %
Max relative error limit: 1.100 %

Test Successful: max relative error is not greater than limit: 1.008% <= 1.100%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_16-49-28/channel_flow_non_newtonian_output
```

![Channel Flow Non-Newtonian Benchmark (Velocity)](images/channel_flow_non_newtonian_velocity_60.png)
##### fig:channel-flow-non-newtonian-velocity
*Benchmark result for velocity using EarthBox version 0.1.0*

![Channel Flow Non-Newtonian Benchmark (Viscosity)](images/channel_flow_non_newtonian_viscosity_60.png)
##### fig:channel-flow-non-newtonian-viscosity
*Benchmark result for effective viscosity using EarthBox version 0.1.0*