# Time Step Initialization

The time step parameters `timestep_viscoelastic` and `ntimestep_max` are initialized based on 
(1) the key word arguments provided by the user to the `run_time_steps` function, (2) whether 
extensional boundary condition have been enabled and (3) whether velocity stepping has been 
enabled.

## Scenario 1: User provided values

The viscoelastic time step (`timestep_viscoelastic`) and number of time step 
(`ntimestep_max`) are based on user provided values via the `run_time_steps` function or loaded from 
the yaml model input file.

## Scenario 2: Symmetric extensional boundary conditions without velocity stepping

The viscoelastic time step (`timestep_viscoelastic`), number of time step (`ntimestep_max`) are 
calculated for models with symmetric extensional boundary conditions without velocity stepping. 

Required conditions for this scenario:
- Time stepping parameters have not been previously set as indicated by `timestep_viscoelastic` and 
  `ntimestep_max` being equal to a default value of NaN.
- Symmetric extensional boundary conditions have been enabled.
- The model duration parameter (`model_duration_myr`) has been set.
- Velocity stepping has not been enabled as indicated by the `iuse_velocity_step` parameter being set to 0.

If the above conditions are met, the viscoelastic time step (`timestep_viscoelastic`) is calculated
using the full extension velocity (`full_extension_velocity`), the maximum allowed cell displacement 
fraction (`marker_cell_displ_max`) and the minimum calculated cell spacing:

###### eq:viscoelastic-timestep-calculated
```math
    \Delta t_{ve} = 
        \min \left( 
            \frac{ \min \left( \Delta x_{min}, \Delta y_{min} \right) f_{cell} }{ 0.5 v_{ext} }, \Delta t_{limit}
        \right)
```

where ``\Delta t_{ve}`` is the viscoelastic time step, ``\Delta x_{min}`` and ``\Delta y_{min}`` are the minimum 
cell spacings x- and y-directions, respectively, ``f_{cell}`` is the maximum allowed cell displacement fraction, 
``v_{ext}`` is the full extension velocity, and ``\Delta t_{limit}`` is the maximum allowed 
viscoelastic time step (50,000 years).

The number of time steps (`ntimestep_max`) is calculated using the viscoelastic time step 
(`timestep_viscoelastic`) and the model duration (`model_duration_myr`) taking a buffer time of 3 Myr 
into account for the adaptive time stepping:

###### eq:number-of-timesteps-calculated
```math
n_{ts,max} = \frac{ \Delta t_{model} + \Delta t_{buffer} }{ \Delta t_{ve} }
```

where ``n_{ts,max}`` is the maximum number of time steps, ``\Delta t_{model}`` is the model duration, 
``\Delta t_{ve}`` is the viscoelastic time step, and ``\Delta t_{buffer}`` is the buffer time.

!!! warning "adaptive time stepping and actual model duration"
    If the calculated magnitude of velocity exceeds the half-extension rate, the viscoelastic time step 
    may be reduced more than estimated using [Eq. ](@ref eq:viscoelastic-timestep-calculated). Therefore,
    [Eq.](@ref eq:number-of-timesteps-calculated) may be too small for achieving the desired
    model duration. To mitigate this issue set the `model_duration_myr` parameter to a larger
    value.  

## Scenario 3: Symmetric extensional boundary conditions with velocity stepping

The viscoelastic time steps (`timestep_viscoelastic`, `timestep_viscoelastic_step1` and 
`timestep_viscoelastic_step2`), number of time steps (`ntimestep_max`), velocity stepping factors 
(`velocity_step_factor` and `velocity_second_step_factor`) and time step adjustment factors 
(`timestep_adjustment_factor` and `timestep_second_adjustment_factor`) are calculated for models with 
symmetric extensional boundary conditions with velocity stepping. 

Required conditions for this scenario:
- Symmetric extensional boundary conditions have been enabled.
- Velocity stepping factors have not been previously set as indicated by the `velocity_step_factor` equal 
  to the default value of NaN.
- The time stepping parameters have not been previously set as indicating by `timestep_viscoelastic` 
  and `ntimestep_max` being equal to a default value of NaN.
- The model duration parameter (`model_duration_myr`) has been set.
- Velocity stepping has been enabled by setting the `iuse_velocity_step` parameter to 1 
  (see [`BoundaryConditions.VelocityStep.initialize!`](@ref EarthBox.BoundaryConditions.VelocityStep.initialize!)
  for more details).
- The first velocity step time (`velocity_first_step_time`) has been set using the 
  [`BoundaryConditions.VelocityStep.initialize!`](@ref EarthBox.BoundaryConditions.VelocityStep.initialize!) 
  function or loaded from the model input file. The second velocity step time (`velocity_second_step_time`)
  is optional and only required if a second velocity step is desired.  
- Full extensions velocities for all phases (`full_velocity_extension`, `full_velocity_extension_step1` 
  have been set using the [`BoundaryConditions.Velocity.initialize!`](@ref EarthBox.BoundaryConditions.Velocity.initialize!) 
  function or loaded from the model input file. The second full extension velocity 
  (`full_velocity_extension_step2`) is optional and only required if a second velocity step is desired.  

The calculation procedure for this scenario is similar to Scenario 2 but modified to account for 
the multiple time steps, velocity stepping factors and time step adjustment factors associated with 
velocity stepping.