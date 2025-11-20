# Viscoelastic Extension

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

# Symmetric

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:viscoelastic_extension);
```

*From test_results.yml*
```
relative error limit %: 3.0
max relative error percentage: 0.2764
message: "Test Successful: max relative error is not greater than limit: 0.276% <= 3.000%"
result: "Success"
```

![Viscoelastic Extension Benchmark](images/viscoelastic_extension_topo_40.png)
##### fig:viscoelastic-extension
*Benchmark result using EarthBox version 0.1.0*


# Asymmetric

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:viscoelastic_extension_asymmetric);
```

*From test_results.yml*
```
relative error limit %: 3.0
max relative error percentage: 0.2243
message: "Test Successful: max relative error is not greater than limit: 0.224% <= 3.000%"
result: "Success"
```

![Viscoelastic Extension Benchmark (Asymmetric)](images/viscoelastic_extension_asymmetric_topo_40.png)
##### fig:viscoelastic-extension-asymmetric
*Benchmark result using EarthBox version 0.1.0*


# Inflow and Outflow Along Side Boundaries

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:viscoelastic_extension_inflow_and_outflow_along_sides);
```

*From test_results.yml*
```
relative error limit %: 3.0
max relative error percentage: 0.5077
message: "Test Successful: max relative error is not greater than limit: 0.508% <= 3.000%"
result: "Success"
```

![Viscoelastic Extension Inflow/Outflow Along Sides Benchmark](images/viscoelastic_extension_inflow_and_outflow_along_sides_topo_40.png)
##### fig:viscoelastic-extension-inflow-outflow
*Benchmark result using EarthBox version 0.1.0*


## Depth-dependent

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:viscoelastic_extension_depth_dependent);
```

*From test_results.yml*
```
relative error limit %: 3.0
max relative error percentage: 0.2744
message: "Test Successful: max relative error is not greater than limit: 0.274% <= 3.000%"
result: "Success"
```

![Viscoelastic Extension Depth-dependent Benchmark](images/viscoelastic_extension_depth_dependent_topo_40.png)
##### fig:viscoelastic-extension-depth-dependent
*Benchmark result using EarthBox version 0.1.0*
