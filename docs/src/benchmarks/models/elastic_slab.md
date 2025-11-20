# Elastic Slab

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:elastic_slab);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Working on time step 160

Max difference (num-ana): 1503.78021
Max relative error: 0.200 %
Max relative error limit: 0.300 %

Test Successful: max relative error is not greater than limit: 0.200% <= 0.300%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_17-23-35/elastic_slab_output
```

![Elastic Slab Benchmark](images/elastic_slab_distance_from_initial_160.png)
##### fig:elastic-slab
*Benchmark result using EarthBox version 0.1.0*