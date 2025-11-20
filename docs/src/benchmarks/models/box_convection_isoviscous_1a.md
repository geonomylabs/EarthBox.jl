# Box Convection Isoviscous 1a

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:box_convection_isoviscous_1a);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Max difference (num-ana): 0.01676                                              
Max relative error: 0.039 %                                                    
Max relative error limit: 0.200 %                                              
                                                                               
Test Successful: max relative error is not greater than limit: 0.039% <= 0.200%
                                                                               
Look at the following path for benchmark plots: ../earthbox_benchmark_results_2025-11-15_17-48-48/box_convection_isoviscous_1a_output
```

![Plasticity Benchmark](images/box_convection_isoviscous_1a.png)
##### fig:box-convection-isoviscous-1a
*Benchmark result using EarthBox version 0.1.0*
