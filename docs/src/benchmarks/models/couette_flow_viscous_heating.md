# Couette Flow with Viscous Heating

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:couette_flow_viscous_heating);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 3.38579
Max relative error: 0.280 %
Max relative error limit: 0.300 %

Test Successful: max relative error is not greater than limit: 0.280% <= 0.300%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_12-52-30/couette_flow_viscous_heating_output
```

![Couette Flow Viscous Heating Benchmark](images/couette_flow_viscous_heating_100.png)
##### fig:couette-flow-viscous-heating
*Benchmark result using EarthBox version 0.1.0*