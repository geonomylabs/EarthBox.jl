# Viscoelastic Stress Buildup

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:viscoelastic_stress_buildup);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 0.00655
Max relative error: 0.076 %
Max relative error limit: 0.100 %

Test Successful: max relative error is not greater than limit: 0.076% <= 0.100%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_17-37-19/viscoelastic_stress_buildup_output
```

![Viscoelastic Stress Buildup Benchmark](images/viscoelastic_stress_buildup_25.png)
##### fig:viscoelastic-stress-buildup
*Benchmark result using EarthBox version 0.1.0*

