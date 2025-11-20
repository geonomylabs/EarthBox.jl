# MarkerTemperature

# Initial Temperature Models

```@eval
using EarthBox
import Markdown
options = Markers.MarkerTemperature.get_options()
liststr = join(["- $(options[id].option_name)" for id in keys(options)], "\n")
Markdown.parse(liststr)
```

!!! tip "quick search for initial temperature models"
    Highlight an initial temperature model name in the list above and use `Ctl-F` + `Enter`
    to navigate to a detailed description in the **Initial Marker Temperature Models** section below.

# Initialization

```@docs
Markers.MarkerTemperature.initialize!
```

# Example: AnalyticalThreeLayer

```julia

# Copy and paste this code into the Julia REPL.
using EarthBox
using CairoMakie
# Initialize EarthBox including model data container.
eb = EarthBoxState(
    xsize=500000.0, ysize=160000.0, 
    xnum=250, ynum=80, 
    dx_marker=500.0, dy_marker=500.0);
# Unpack the model data container.
model = eb.model_manager.model;
# Initialize the staggered grid.
StaggeredGrid.initialize!(model, grid_type=:UniformGrid);
# Initialize marker coordinates.
Markers.MarkerCoordinates.initialize!(model, marker_distribution=:Regular);
# Initialize the required material geometry parameters for the 
# LithosphericExtensionLateralStrongZones material model
MaterialGeometry.StickyAirGeometry.initialize!(model, thick_air=10000.0);
MaterialGeometry.EarthLayering.initialize!(
    model, thick_lith=125000.0, thick_crust=35000.0, thick_upper_crust=22000.0);
# Get EarthBox parameters to programmatically define parameter dictionary keys
const PARAMS = get_eb_parameters();
# Initialize marker temperature
Markers.MarkerTemperature.initialize!(
    model,
    initial_temperature_model=initial_temperature_names.AnalyticalThreeLayer,
    parameters=Dict(
        PARAMS.temperature_top.name             => ConversionFuncs.celsius_to_kelvin(0.0),
        PARAMS.amplitude_perturbation.name      => 0.0,
        PARAMS.width_perturbation.name          => 10_000.0,
        PARAMS.thick_thermal_lithosphere.name   => 125_000.0,
        PARAMS.temperature_base_lith.name       => ConversionFuncs.celsius_to_kelvin(1345.0),
        PARAMS.adiabatic_gradient.name          => 0.4,
        PARAMS.conductivity_upper_crust.name    => 2.25,
        PARAMS.conductivity_lower_crust.name    => 2.0,
        PARAMS.conductivity_mantle.name         => 2.0,
        PARAMS.heat_production_upper_crust.name => 1.8e-6,
        PARAMS.heat_production_lower_crust.name => 0.5e-6,
        PARAMS.heat_production_mantle.name      => 0.0
    )
);
# Get marker coordinate arrays amd convert from meters to kilometers
marker_x_km = model.obj_dict["marker_x"].array./1000.0;
marker_y_km = model.obj_dict["marker_y"].array./1000.0;
# Get marker temperature array and convert to Celsius
marker_TC = ConversionFuncs.kelvin_to_celsius.(model.obj_dict["marker_TK"].array);
# Get size parameters and convert from meters to kilometers
xsize_km = model.obj_dict["xsize"].value/1000.0;
ysize_km = model.obj_dict["ysize"].value/1000.0;
# Decimate the marker arrays for faster plotting
decimation_factor = 2;
indices = 1:decimation_factor:length(marker_TC);
marker_TC = marker_TC[indices];
marker_x_km = marker_x_km[indices];
marker_y_km = marker_y_km[indices];
# Plot using CairoMakie
dpi = 150.0;
fig = Figure(size = (10.0*dpi, 3.0*dpi)); # 10 inches x 3 inches @ 150 dpi
ax = Axis(
    fig[1, 1]; limits = (0.0, xsize_km, 0.0, ysize_km), yreversed = true,
    xticks = collect(0.0:25.0:xsize_km), yticks = collect(0.0:10.0:ysize_km),
    xlabel = "X (km)", ylabel = "Y (km)", 
    title = "AnalyticalThreeLayer Model for Initial Marker Temperature", 
    titlesize = 12, xlabelsize = 12, ylabelsize = 12, 
    xticklabelsize = 12, yticklabelsize = 12
);
plt = scatter!(
    ax, marker_x_km, marker_y_km, color = marker_TC, 
    colormap = :bwr,
    colorrange = (0.0, 1350.0),
    marker=:circle, markersize = 4, rasterize = true
); 
cb = Colorbar(
    fig[1,2], plt; label = "Marker Temperature (Celsius)",
    ticks = collect(0.0:50.0:1350.0),
    ticklabelsize = 10
);
fig

```