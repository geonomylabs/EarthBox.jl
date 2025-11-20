# Batch Benchmark Execution

```@docs
BenchmarksManager.run_benchmarks
```
    
## Example: Running Multiple "Fast" Benchmark Cases

This example shows how multiple benchmarks can be run in batch mode with output
being sent to the present working directory. Check terminal output for a quick look
at test results and for the path to where results were written to file.

```julia

using DataStructures: OrderedDict
using EarthBox

test_dict = OrderedDict(
    benchmark_names.couette_flow_viscous_heating        => true,
    benchmark_names.channel_flow_non_newtonian          => true,
    benchmark_names.channel_flow_variable_conductivity  => true,
    benchmark_names.channel_flow_non_steady_temperature => true,
    benchmark_names.solid_body_rotation                 => true,
    benchmark_names.rayleigh_taylor_instability         => true,
);
BenchmarksManager.run_benchmarks(test_dict);

```

This code will produce an output directory with a time-stamped name like 
`earthbox_benchmark_results_2025-11-07_08-48-04` that contains model output,
benchmark plots and a summary file of test results called `test_results.yml`.


## Example: Running Multiple Benchmark Cases with User Configuration

This example shows how multiple benchmarks can be run in batch mode with control
over model execution, post processing, output base path and the use of the MUMPS
solver.

```julia

using DataStructures: OrderedDict
using EarthBox

benchmark_names = BenchmarksManager.benchmark_names
test_dict = OrderedDict(
    benchmark_names.couette_flow_viscous_heating                          => true,
    benchmark_names.channel_flow_non_newtonian                            => true,
    benchmark_names.channel_flow_variable_conductivity                    => true,
    benchmark_names.channel_flow_non_steady_temperature                   => true,
    benchmark_names.solid_body_rotation                                   => true,
    benchmark_names.rayleigh_taylor_instability                           => true,
    benchmark_names.elastic_slab                                          => true,
    benchmark_names.viscoelastic_stress_buildup                           => true,
    benchmark_names.box_convection_isoviscous_1a                          => true,
    benchmark_names.plasticity_benchmark_kaus10                           => true,
    benchmark_names.viscoelastic_extension                                => false,
    benchmark_names.viscoelastic_extension_asymmetric                     => false,
    benchmark_names.viscoelastic_extension_depth_dependent                => false,
    benchmark_names.viscoelastic_extension_inflow_and_outflow_along_sides => false,
    benchmark_names.viscoelastic_contraction                              => false,
    benchmark_names.viscoelastic_contraction_asymmetric                   => false,
    benchmark_names.simple_sedimentation                                  => false,
    benchmark_names.seafloor_spreading                                    => false,
    benchmark_names.flexure_triangular_hole                               => true
)
# Dictionary used to activate the mumps solver for large models
# For each dictionary key the value is a vector with two elements:
# [1] Bool: Activate/deactivate mumps solver
# [2] Int: Number of processors to use
mumps_solver_dict = Dict{Symbol, Vector{Any}}(
    benchmark_names.flexure_triangular_hole => [true, 8]
)
BenchmarksManager.run_benchmarks(
    test_dict,
    mumps_solver_dict   = mumps_solver_dict,
    run_model           = true,
    run_post_processing = true,
    base_path           = "/mnt/extradrive1",
    make_backup         = false,
    restart_from_backup = false
)

```