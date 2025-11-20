# RockProperties

# Initialization

```@docs
RockProperties.initialize!
```

# Example

```julia
using EarthBox
# Initialize EarthBox including model data container.
eb = EarthBoxState(
    xsize=500000.0, ysize=160000.0, 
    xnum=250, ynum=80, 
    dx_marker=500.0, dy_marker=500.0);
# Unpack the model data container.
model = eb.model_manager.model;
# Initialize a rock property model with the density model from Liao14,
# temperature-dependent thermal conductivity and a rhocp (Density x Heat Capacity) 
# model with temperature dependent heat capacity
RockProperties.initialize!(
    model,
    density_model              = density_model_names.Liao14,
    thermal_conductivity_model = thermal_conductivity_model_names.SekiguchiWaples,
    rhocp_model                = rhocp_model_names.TemperatureDependentWaples,
    iuse_sed_porosity          = 1
);
RockProperties.print_option(model);
```