# Examples

## Example 1: User provided values

```julia
using EarthBox

eb = EarthBoxState(
    xnum = 101, ynum = 101, 
    xsize = 100000.0, ysize = 100000.0, 
    dx_marker = 25.0, dy_marker = 25.0
);

# Additional API calls are required to initialize components of the model

run_time_steps(
    eb,
    make_backup           = true, 
    ntimestep_max         = 4500,
    timestep_viscoelastic = years_to_seconds(50_000.0), 
    timestep_out          = years_to_seconds(500_000.0)
    );
```

## Example 2: Symmetric extensional boundary conditions without velocity stepping

This example demonstrates how to calculate derived parameters for testing 
purposes without the need to initialize a full model with additional API calls. 
Replace TimeLoop.initialize! with a call to `run_time_steps` after a model is 
fully initialized with additional API calls to derive parameters and run a full
model. Note that the `run_time_steps` function performs time loop initialization
internally so the `TimeLoop.initialize!` call is not required.

```julia
using EarthBox

eb = EarthBoxState(
    xnum = 101, ynum = 101, 
    xsize = 100000.0, ysize = 100000.0, 
    dx_marker = 100.0, dy_marker = 100.0);

model = eb.model_manager.model;

BoundaryConditions.initialize!(
    model, 
    model_type = "LithosphericExtensionFixedBoundaries");

BoundaryConditions.Velocity.initialize!(
    model, full_velocity_extension = ConversionFuncs.cm_yr_to_m_s(4.0));

TimeLoop.initialize!(
    model,
    model_duration_myr = 30.0,
    timestep_out = ConversionFuncs.years_to_seconds(500_000.0));
```

## Example 3: Symmetric extensional boundary conditions with velocity stepping

This example demonstrates how to calculate derived parameters for testing 
purposes without the need to initialize a full model with additional API calls. 
Replace TimeLoop.initialize! with a call to `run_time_steps` after a model is 
fully initialized with additional API calls to derive parameters and run a full
model. Note that the `run_time_steps` function performs time loop initialization
internally so the `TimeLoop.initialize!` call is not required.

```julia
using EarthBox

eb = EarthBoxState(
    xnum = 101, ynum = 101, 
    xsize = 100000.0, ysize = 100000.0, 
    dx_marker = 100.0, dy_marker = 100.0);

model = eb.model_manager.model;

BoundaryConditions.initialize!(
    model, 
    model_type = "LithosphericExtensionFixedBoundaries");

BoundaryConditions.Velocity.initialize!(
    model, 
    full_velocity_extension = ConversionFuncs.cm_yr_to_m_s(1.0), 
    full_velocity_extension_step1 = ConversionFuncs.cm_yr_to_m_s(1.0), 
    full_velocity_extension_step2 = ConversionFuncs.cm_yr_to_m_s(4.0));

BoundaryConditions.VelocityStep.initialize!(
    model, 
    iuse_velocity_step = 1, 
    velocity_step_time = 10.0, 
    velocity_second_step_time = 20.0);

TimeLoop.initialize!(
    model,
    model_duration_myr = 16.0,
    timestep_out = ConversionFuncs.years_to_seconds(500_000.0));
```