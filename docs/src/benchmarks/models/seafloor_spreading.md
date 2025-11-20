# Seafloor Spreading

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:seafloor_spreading);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
Average off-axis thickness of oceanic crust (Numerical): 7405.106679485966
Average thickness of oceanic crust (Target): 7000.0                            
                                                                               
Max difference (num-ana): 405.10668                                            
Max relative error: 5.787 %                                                    
Max relative error limit: 6.000 %                                              
                                                                               
Test Successful: max relative error is not greater than limit: 5.787% <= 6.000%
```

![Seafloor Spreading Benchmark (Topography)](images/seafloor_spreading_topo_60.png)
##### fig:seafloor-spreading-topography
*Benchmark result for topography using EarthBox version 0.1.0*

![Seafloor Spreading Benchmark (Magmatic Crustal Thickness)](images/seafloor_spreading_xth_60.png)
##### fig:seafloor-spreading-xth
*Benchmark result for magmatic crustal thickness using EarthBox version 0.1.0*


