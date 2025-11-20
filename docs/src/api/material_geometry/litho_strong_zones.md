# LithoStrongZones

The `LithoStrongZones` geometry is composed of a central upper weak zone in the lithosphere
surrounded by strong lithospheric zones.

# Initialization

```@docs
ModelManager.MaterialGeometry.LithoStrongZones.initialize!
```

# Initialization Scenarios

Initialization involves defining the left and right boundaries of the central 
lithospheric zone, `x_left_strong` and `x_right_strong`, respectively. Also, the 
flag `iuse_strong_zones` argument must be set to 1 to enable strong zones.

Initialization of these parameters depends on whether T-type refinement 
parameters `xo_highres` and `xf_highres` were previously defined during model 
initialization with [`EarthBoxState`](@ref EarthBoxState) or during staggered 
grid initialization with [`StaggeredGrid.initialize!`](@ref EarthBox.StaggeredGrid.initialize!).

If T-type refinement parameters `xo_highres` and `xf_highres` were previously 
defined then `x_left_strong` and `x_right_strong` of the central lithospheric 
zone will be set equal to `xo_highres` and `xf_highres`, respectively.

If T-type refinement parameters `xo_highres` and `xf_highres` were not previously 
defined then `x_left_strong` and `x_right_strong` will be set to equal to the 
values provided as arguments to `initialize!`.

# Examples

## Example 1: Direct input of strong zone boundaries

```julia
eb = EarthBoxState(xnum = 100, ynum = 100, xsize = 500000.0, ysize = 120000.0)
model = eb.model_manager.model # ModelData container
EarthBox.MaterialGeometry.LithoStrongZones.initialize!(
    model, 
    iuse_strong_zones = 1,
    x_left_strong = 200000.0,
    x_right_strong = 300000.0
)
```

## Example 2: Using T-type refinement parameters to set strong zone boundaries

```julia
ttype_refinement_parameters = Dict(
    "xo_highres" => 200000.0,
    "xf_highres" => 300000.0,
    "yf_highres" => 50000.0,
    "dx_highres" => 500.0,
    "dy_highres" => 500.0,
    "dx_lowres" => 2000.0,
    "dy_lowres" => 2000.0
)
eb = EarthBoxState(
    xsize = 500000.0, ysize = 120000.0, 
    ttype_refinement_parameters = ttype_refinement_parameters
)
model = eb.model_manager.model
EarthBox.MaterialGeometry.LithoStrongZones.initialize!(model, iuse_strong_zones = 1)
```