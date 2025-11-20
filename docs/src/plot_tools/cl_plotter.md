# Command-line Plotter

The EarthBox command-line plotter allows the user to build customized functions 
for plotting EarthBox output files. The standard way to use the EarthBox plotter 
involves creating a script called `Plot.jl` in a directory containing model input 
files and scripts. The `Plot.jl` script should contain a module called `Plot` with 
plotting functions, model output paths, material definition (for marker and yield plots) 
and a function that executes the plotter function `run_cl_plotter`. Paths and options
are specified via command line with an option to directly pass the root path of the
model output directory as a function argument to eliminate the need to write the 
full path via command line.

For more detailed information about available plotting functions that can be used 
with the command-line plotter see the following:
- [Plotting Scalar Arrays](@ref)
- [Plotting Marker Swarms](@ref)
- [Plotting Velocity Vectors](@ref)
- [Plotting Yield Strength](@ref)
- [Plotting Stokes Convergence](@ref)

```@docs
ModelPlots2DManager.run_cl_plotter
```