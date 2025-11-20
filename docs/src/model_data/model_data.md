# ModelData
```@docs

ModelDataContainer.ModelData

```

# Accessing ModelData Objects

Parameter and array objects can be access from the [ModelData](@ref) structure
via the `obj_dict` field.

```julia

using EarthBox

eb = EarthBoxState(xnum=100, ynum=100, xsize=100000, ysize=100000.0, dx_marker=100.0, dy_marker=100.0);

model = eb.model_manager.model;
xnum = model.obj_dict["xnum"];
marker_x = model.obj_dict["marker_x"];
marker_y = model.obj_dict["marker_y"];

x = marker_x.array[50]
array_units = marker_x.units
array_description = marker_x.description

value = xnum.value
units = xnum.units
parameter_description = xnum.description

```

!!! note "nested dot access"
    Model objects can also be accessed via nested dot structures. For example, the marker x-coordinates
    can be accessed via `model.markers.arrays.location.marker_x.array`. The nested dot access patterns
    are defined for each objects in a model data collection. See 
    [MarkerContainer.ArrayCollection.Location][@ref ModelDataContainer.MarkerContainer.ArrayCollection.Location]
    for an example of nested dot access definitions.
