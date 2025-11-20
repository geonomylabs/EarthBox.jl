# Channel Flow Non-steady Temperature

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:channel_flow_non_steady_temperature);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 4.86148
Max relative error: 0.658 %
Max relative error limit: 0.700 %

Test Successful: max relative error is not greater than limit: 0.658% <= 0.700%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_16-39-32/channel_flow_non_steady_temperature_output
```

![Channel Flow Non-steady Temperature](images/channel_flow_non_steady_temperature_temperature_240.png)
##### fig:channel-flow-non-steady-temperature
*Benchmark result using EarthBox version 0.1.0*

