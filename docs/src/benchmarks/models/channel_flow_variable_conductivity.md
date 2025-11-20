# Channel Flow Variable Conductivity

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:channel_flow_variable_conductivity);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 1.83562
Max relative error: 2.105 %
Max relative error limit: 2.200 %

Test Successful: max relative error is not greater than limit: 2.105% <= 2.200%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_14-53-27/channel_flow_variable_conductivity_output
```

![Couette Flow Viscous Heating Benchmark](images/channel_flow_variable_conductivity_temperature_100.png)
##### fig:channel-flow--variable-conductivity
*Benchmark result using EarthBox version 0.1.0*

