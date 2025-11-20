# Flexure Triangular Hole

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:flexure_triangular_hole);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
relative error limit %: 15.0
max relative error percentage: 13.2312
message: "Test Successful: max relative error is not greater than limit: 13.231% <= 15.000%"
result: "Success"
```

![Flexure Triangular Hole Benchmark](images/flexure_triangular_hole_topo_60.png)
##### fig:flexure-triangular-hole
*Benchmark result using EarthBox version 0.1.0*

