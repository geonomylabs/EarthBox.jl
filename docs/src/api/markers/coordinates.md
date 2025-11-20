# MarkerCoordinates

# Initialization

```@docs
Markers.MarkerCoordinates.initialize!
```

# Accessing Option Names

Available marker distribution options can be accessed programmatically:
- `Markers.MarkerCoordinates.option_names.Regular`
- `Markers.MarkerCoordinates.option_names.Randomized`


# Example: Initializing a model with a randomized marker distribution

```julia
using EarthBox
using CairoMakie

eb = EarthBoxState(
    xnum = 100, ynum = 100, 
    xsize = 100000.0, ysize = 100000.0, 
    dx_marker = 1000.0, dy_marker = 1000.0);

model = eb.model_manager.model;

option_names = Markers.MarkerCoordinates.option_names;
Markers.MarkerCoordinates.initialize!(model, marker_distribution=option_names.Randomized);

marker_x = model.markers.arrays.location.marker_x.array;
marker_y = model.markers.arrays.location.marker_y.array;

fig, ax, plt = scatter(
    marker_x ./ 1000.0, marker_y ./ 1000.0; 
    color = :dodgerblue, markersize = 2, 
    axis = (; aspect = 1, xlabel = "X (km)", ylabel = "Y (km)", limits=(0.0, 100.0, 0.0, 100.0)),
    );
fig

```