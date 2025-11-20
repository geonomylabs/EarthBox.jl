# Solid Body Rotation

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:solid_body_rotation);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 0.07100
Max relative error: 0.010 %
Max relative error limit: 2.000 %

Test Successful: max relative error is not greater than limit: 0.010% <= 2.000%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_17-14-19/solid_body_rotation_output
```

![Solid Body Rotation Benchmark](images/solid_body_rotation_temperature_271.png)
##### fig:solid-body-rotation
*Benchmark result using EarthBox version 0.1.0*