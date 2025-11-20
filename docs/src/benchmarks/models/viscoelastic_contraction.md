# Viscoelastic Contraction

See descriptions at [EarthBox.BenchmarksManager.run_benchmark](@ref).

# Symmetric

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:viscoelastic_contraction);
```

*From test_results.yml*
```
relative error limit %: 3.0
max relative error percentage: 2.0714
message: "Test Successful: max relative error is not greater than limit: 2.071% <= 3.000%"
result: "Success"
```
![Viscoelastic Contraction](images/viscoelastic_contraction_topo_8.png)
##### fig:viscoelastic-contraction
*Benchmark result using EarthBox version 0.1.0*


## Asymmetric

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:viscoelastic_contraction_asymmetric);
```

*From test_results.yml*
```
relative error limit %: 3.0
max relative error percentage: 0.7165
message: "Test Successful: max relative error is not greater than limit: 0.717% <= 3.000%"
result: "Success"
```

![Viscoelastic Contraction (Asymmetric)](images/viscoelastic_contraction_asymmetric_topo_8.png)
##### fig:viscoelastic-contraction-asymmetric
*Benchmark result using EarthBox version 0.1.0*
