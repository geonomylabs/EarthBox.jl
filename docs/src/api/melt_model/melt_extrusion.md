# MeltExtrusion

# Initialization

Extruded material is defined using material types called 
"ExtrudedGabbroicMagma" (extruded lava) and "SolidifiedBasalt"
(solidified lava). These two material types must be defined in the 
material model for the extrusion model to work. See [Material Types](@ref)
for more information about material types. 

The volume of extruded volcanic material is calculated as the product of
the extrusion volume factor and the total incremental volume of magma 
produced within the mantle in a given migration domain.

The flow distance of individual lava flows is controlled by the residual lava
thickness, total volume of the flow and surface relief. Lava flows will 
advance only after the residual lava thickness has been reached and a 
downward slope is present.

The characteristic volume of lava per flow is equal to the product of
the characteristic flow length and the residual lava thickness, which 
have different values for subaerial and submarine flows. The total
number of lava flows per time step is calculated as the total volume of
magma produced in a particular migration domain divided by the characteristic
volume of lava per flow. The flow distance of a given lava flow is
directly proportional to the characteristic flow length and inversely 
proportional to the residual lava thickness.

The location of extrusion at the surface is controlled by the extrusion
window. The center of the extrusion window is determined by 
x-coordinate of the shallowest partially molten mantle marker. If 
`width_eruption_domain_fixed` is zero then the width of the extrusion 
window is defined by the minimum and maximum x-coordinates of the 
partially molten or molten gabbroic material. With this model, extrusion 
is spatially related to gabbroic magma chambers and gabbroic magma mush 
zones (gabbro glacier). If `width_eruption_domain_fixed` is greater than
zero, the width of the extrusion window is equal to the value of
`width_eruption_domain_fixed`.

Two models are available for eruption location of a laval flow at the 
surface within the extrusion window:
- Probabilistic:
    - The eruption location is chosen either randomly 
       (`iuse_random_eruption_location = 1`) or from a normal distribution 
       (`iuse_normal_eruption_location = 1`) within the extrusion window.
- Static (`iuse_random_eruption_location = 0`): 
    - With this option the eruption location is fixed at the midpoint of the 
       extrusion window.

The parameters `porosity_initial_lava_flow` and `decay_depth_lava_flow` are
used to decompact extruded gabbroic magma to account for vesicles. Note
that compaction during burial is not controlled by these parameters and
is instead controlled by the material type.

```@docs
MeltModel.Extrusion.initialize!
```