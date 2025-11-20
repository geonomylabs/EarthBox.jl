# Plasticity Benchmark Kaus10

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:plasticity_benchmark_kaus10);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 1.00125                                              
Max relative error: 1.780 %                                                    
Max relative error limit: 4.000 %                                              
                                                                               
Test Successful: max relative error is not greater than limit: 1.780% <= 4.000%
                                                                               
Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_17-48-48/plasticity_benchmark_kaus10_output
```

![Plasticity Benchmark](images/plasticity_benchmark_kaus10.png)
##### fig:plasticity-benchmark-kaus10
*Benchmark result using EarthBox version 0.1.0*

