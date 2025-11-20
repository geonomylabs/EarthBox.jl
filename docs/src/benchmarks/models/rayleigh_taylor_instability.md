# Rayleigh Taylor Instability

See description at [EarthBox.BenchmarksManager.run_benchmark](@ref).

This benchmark can be run with the following code:

```julia
using EarthBox
BenchmarksManager.run_benchmark(:rayleigh_taylor_instability);
```

Output will be sent to a directory in the present working directory. Output will
also be sent to the terminal indicating the result of the benchmark and
where to find output files.

```
>> Rayleigh-Taylor Instability Benchmark Input and Output
   layer1_thickness: 1500.0000
   layer2_thickness: 1500.0000
   density_layer1: 3000.0000
   density_layer2: 2900.0000
   viscosity_layer1: 10000000000000.0000
   viscosity_layer2: 10000000000000000000.0000
   amplitude_initial: 10.0000
   wavelength: 4000.0000
   growth_factor_k_numerical: 0.3137
   velocity_y_wave_numerical(cm/yr): -7.2836E-04
   growth_factor_k_analytical: 0.3190
   velocity_y_wave_analytical(cm/yr): -7.4077E-04
   theta1: 2.3562
   b1: 0.5000
   b2: 0.2000
   b1*K+b2 (numerical): 0.3568
   b1*K+b2 (analytical): 0.3595


Max difference (num-ana): 0.00000
Max relative error: 1.675 %
Max relative error limit: 3.500 %

Test Successful: max relative error is not greater than limit: 1.675% <= 3.500%

Look at the following path for benchmark plots: .../earthbox_benchmark_results_2025-11-15_17-18-47/rayleigh_taylor_instability_output

```

