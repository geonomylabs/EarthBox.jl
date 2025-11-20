# MarkerMaterials

The `MarkerMaterials` module is used to initialize the material model used in the
EarthBox simulation.

# Material Models

```@eval
using EarthBox
import Markdown
options = Markers.MarkerMaterials.get_options()
liststr = join(["- $(options[id].option_name)" for id in keys(options)], "\n")
Markdown.parse(liststr)
```

!!! tip "quick search for material model"
    Highlight a material model name in the list above and use `Ctl-F` + `Enter`
    to navigate to a detailed description in the **Material Models** section below.

# Initialization

```@docs
Markers.MarkerMaterials.initialize!
```

# Example

This example illustrates how to initialize marker material definitions for a case
using the `LithosphericExtensionLateralStrongZones` material model with an upper 
and lower continental crust, layered continental mantle lithosphere, and lateral 
zones with stronger frictional plastic properties. Note that this material model
requires certain material domains and material geometry input parameters. See 
description provided for `LithosphericExtensionLateralStrongZones` in the 
`Material Models` section of [MarkerMaterials.initialize!](@ref Markers.MarkerMaterials.initialize!) 
for more details. Also see [Material Domains](@ref) and [Material Types](@ref)
for more information about material domains and types, respectively.

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
# Define the material collection with material names and associated properties
collection =  MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1;
# You can inspect the material collection being used by opening the collection
# file associated with the following path. You can also make a local copy of this
# file and create a custom material collection. For this custom case the you will
# need to update the custom_path variable to your custom location.
collection_path = collection.path
# Unpack material names from the collection, names for material types and names
# for material domains.
mat_names = collection.materials;
types = MaterialTypesRegistry();
domains = MaterialDomainsRegistry();
# Define the material input taking into account the required material domains
# for the LithosphericExtensionLateralStrongZones material model that will be
# initialize during material initialization.
materials_input_dict = MaterialsDictType(
    # Sticky Domain
    Int16(1) => MaterialDictType(
        "mat_name" => mat_names.sticky_air,
        "mat_type" => types.sticky_air,
        "mat_domain" => domains.atmosphere,
        "red_fraction" => 255/255,
        "green_fraction" => 255/255,
        "blue_fraction" => 255/255,
    ),
    Int16(2) => MaterialDictType(
        "mat_name" => mat_names.sticky_water,
        "mat_type" => types.sticky_water,
        "mat_domain" => domains.ocean,
        "red_fraction" => 0/255,
        "green_fraction" => 255/255,
        "blue_fraction" => 255/255,
    ),
    # Felsic Continental Crust
    Int16(3) => MaterialDictType(
        "mat_name" => mat_names.felsic_continental_crust,
        "mat_type" => types.felsic_continental_crust_fertile,
        "mat_domain" => domains.upper_continental_crust,
        "red_fraction" => 255/255,
        "green_fraction" => 153/255,
        "blue_fraction" => 153/255,
    ),
    Int16(4) => MaterialDictType(
        "mat_name" => mat_names.felsic_continental_crust_strong_zone,
        "mat_type" => types.felsic_continental_crust_fertile,
        "mat_domain" => domains.upper_continental_crust_strong_zone,
        "red_fraction" => 255/255,
        "green_fraction" => 163/255,
        "blue_fraction" => 163/255,
    ),
    # Mafic Continental Crust
    Int16(5) => MaterialDictType(
        "mat_name" => mat_names.mafic_continental_crust,
        "mat_type" => types.mafic_continental_crust_fertile,
        "mat_domain" => domains.lower_continental_crust,
        "red_fraction" => 255/255,
        "green_fraction" => 200/255,
        "blue_fraction" => 200/255,
    ),
    Int16(6) => MaterialDictType(
        "mat_name" => mat_names.mafic_continental_crust_strong_zone,
        "mat_type" => types.mafic_continental_crust_fertile,
        "mat_domain" => domains.lower_continental_crust_strong_zone,
        "red_fraction" => 255/255,
        "green_fraction" => 210/255,
        "blue_fraction" => 210/255,
    ),
    # Continental Mantle Lithosphere
    Int16(7) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_continental_lithosphere,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.upper_mantle_lithosphere,
        "red_fraction" => 0.0,
        "green_fraction" => 153/255,
        "blue_fraction" => 153/255,
    ),
    Int16(8) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_continental_lithosphere,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.middle_mantle_lithosphere,
        "red_fraction" => 0.0,
        "green_fraction" => 153/255,
        "blue_fraction" => (153+25)/255,
    ),
    Int16(9) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_continental_lithosphere,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.lower_mantle_lithosphere,
        "red_fraction" => 0.0,
        "green_fraction" => 153/255,
        "blue_fraction" => 153/255,
    ),
    Int16(10) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_continental_lithosphere_strong_zone,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.lithospheric_mantle_strong_zone,
        "red_fraction" => 0.0,
        "green_fraction" => 163/255,
        "blue_fraction" => 163/255,
    ),
    # Asthenosphere
    Int16(11) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.asthenosphere,
        "red_fraction" => 0.0,
        "green_fraction" => 200/255,
        "blue_fraction" => 153/255,
    )
);
# Initialize the required material geometry parameters for the 
# LithosphericExtensionLateralStrongZones material model
MaterialGeometry.StickyAirGeometry.initialize!(model, thick_air=10000.0);
MaterialGeometry.EarthLayering.initialize!(
    model, thick_lith=125000.0, thick_crust=35000.0, thick_upper_crust=22000.0);
MaterialGeometry.LithoStrongZones.initialize!(
    model, iuse_strong_zones=1, x_left_strong=150000.0, x_right_strong=350000.0);
# Initialize the material model
Markers.MarkerMaterials.initialize!(
    model, 
    material_model       = material_model_names.LithosphericExtensionLateralStrongZones,
    paths                = Dict("materials_library_file" => collection_path),
    materials_input_dict = materials_input_dict,
    viscosity_min        = 1e18,
    viscosity_max        = 1e25 
);
# Get the marker arrays for material ID's and x- and y-coordinates and convert
# from meters to kilometers
marker_matid = model.obj_dict["marker_matid"].array;
marker_x_km = model.obj_dict["marker_x"].array./1000.0;
marker_y_km = model.obj_dict["marker_y"].array./1000.0;
# Get size parameters and convert from meters to kilometers
xsize_km = model.obj_dict["xsize"].value/1000.0;
ysize_km = model.obj_dict["ysize"].value/1000.0;
# Decimate the marker arrays for faster plotting
decimation_factor = 2;
indices = 1:decimation_factor:length(marker_matid);
marker_matid = marker_matid[indices];
marker_x_km = marker_x_km[indices];
marker_y_km = marker_y_km[indices];
# Plot using CairoMakie
dpi = 150.0;
fig = Figure(size = (10.0*dpi, 3.0*dpi)); # 10 inches x 3 inches @ 150 dpi
ax = Axis(
    fig[1, 1]; limits = (0.0, xsize_km, 0.0, ysize_km), yreversed = true,
    xticks = collect(0.0:25.0:xsize_km), yticks = collect(0.0:10.0:ysize_km),
    xlabel = "X (km)", ylabel = "Y (km)", 
    title = "Earth Layering with Lateral Strong Zones", 
    titlesize = 12, xlabelsize = 12, ylabelsize = 12, 
    xticklabelsize = 12, yticklabelsize = 12
);
plt = scatter!(
    ax, marker_x_km, marker_y_km, color = marker_matid, 
    colormap = model.materials.colors.color_map,
    colorrange = model.materials.colors.colorrange,
    marker=:circle, markersize = 4, rasterize = true
);
cb = Colorbar(
    fig[1,2], plt; label = "Marker Materials",
    ticks = (model.materials.colors.ticks, model.materials.colors.labels),
    ticklabelsize = 10
);
fig

```
