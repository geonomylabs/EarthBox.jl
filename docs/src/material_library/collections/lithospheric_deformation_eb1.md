# Material Collection: `lithospheric_deformation_eb1`

# Material List
- [StickyAir EB1](@ref)
- [StickyWater EB1](@ref)
- [ClasticSediment EB1](@ref)
- [Salt EB1](@ref)
- [FelsicContinentalCrust EB1](@ref)
- [FelsicContinentalCrustWeak EB1](@ref)
- [FelsicContinentalCrustStrongZone EB1](@ref)
- [MaficContinentalCrust EB1](@ref)
- [MaficContinentalCrustStrongZone EB1](@ref)
- [UltramaficContinentalLithosphere EB1](@ref)
- [UltramaficContinentalLithosphereWeak EB1](@ref)
- [UltramaficContinentalLithosphereStrongZone EB1](@ref)
- [UltramaficAsthenosphereDryFertile EB1](@ref)
- [UltramaficAsthenosphereWetFertile EB1](@ref)
- [OceanicGabbroicCrust EB1](@ref)
- [GabbroicMagma EB1](@ref)
- [LayeredGabbroicCrust EB1](@ref)
- [LayeredGabbroicMagma EB1](@ref)
- [SerpentinizedPeridotite EB1](@ref)

## Material Summary

**Collection File Path**
```
src/material_library/registries/lithospheric_deformation/lithospheric_deformation_eb1.yaml
```

- Mantle flow laws from HK03.
- Base case crustal flow strength from GT95 and R06 (no scaling factor)
- Moderate plastic yield strength (20/10; 0/0.1)
- Moderate radiogenic heat production (1.8e-6/0.5e-6)
  
### Upper Crustal Rheology (Moderate Strength)

Wet quartzite dislocation creep flow law parameters from Gleason and Tullis (1995). 
No scaling factor is applied to pre-exponential factor.

### Lower Crustal Rheology
Wet anorthite dislocation creep flow law parameters from Rybacki et al. (2006) using a water 
fugacity of 500 MPa based on fit to Figure 8b at 30 km. No scaling factor is applied to 
pre-exponential factor.

### Mantle Rheology
Dislocation and diffusion creep flow law parameters from Hirth and Kohlstedt (2003) for 
dry and wet olivine aggregates. A grain size of 6 mm is assumed for the olivine diffusion creep flow 
law. A water content of 500 ppm H/Si is assumed for the wet olivine aggregates.

Peierls creep law parameters are inferred from Katayama and Karato (2008). The current 
pre-exponential factor for Peierls creep was adjusted from 6.31e+7 (wet case)
to 2.5e+7 approximate dry case and to ensure consistency with the dry dislocation 
creep law. A visual comparison to Figure 9b of Katayama and Karato (2008) was used to
test the adjusted pre-exponential factor at strain rate 1e-5 (~1000 MPa) and 1e-10 (~100 MPa).
Note that the activation energy and volume for Peierls creep comes from the
current dislocation creep law, which requires adjustment to the Peierls 
pre-exponential term if different from Katayama and karato (2008).

### Frictional Plastic Model

Cohesion: 20 to 10 MPa from plastic strain of 0.0 to 1.0.

Friction angles: 30 to 7 degrees from plastic strain of 0.0 to 1.0. The 7 degree minimum
friction angle is from the model presented by Bos and Spiers (2002) based on laboratory 
experiments (see Figure 6 for the case with a sliding velocity of 1e-3 micrometers/second).

### Radiogenic Heat Production
- upper crust: 1.8 microW/m^3
- lower crust: 0.5 microW/m^3

## Material Tables

### StickyAir EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "StickyAir"
    )
Markdown.parse(table_str)
```

### StickyWater EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "StickyWater"
    )
Markdown.parse(table_str)
```

### ClasticSediment EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "ClasticSediment"
    )
Markdown.parse(table_str)
```

### Salt EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "Salt"
    )
Markdown.parse(table_str)
```

### FelsicContinentalCrust EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "FelsicContinentalCrust"
    )
Markdown.parse(table_str)
```

### FelsicContinentalCrustWeak EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "FelsicContinentalCrustWeak"
    )
Markdown.parse(table_str)
```

### FelsicContinentalCrustStrongZone EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "FelsicContinentalCrustStrongZone"
    )
Markdown.parse(table_str)
```

### MaficContinentalCrust EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "MaficContinentalCrust"
    )
Markdown.parse(table_str)
```

### MaficContinentalCrustStrongZone EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "MaficContinentalCrustStrongZone"
    )
Markdown.parse(table_str)
```

### UltramaficContinentalLithosphere EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "UltramaficContinentalLithosphere"
    )
Markdown.parse(table_str)
```

### UltramaficContinentalLithosphereWeak EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "UltramaficContinentalLithosphereWeak"
    )
Markdown.parse(table_str)
```

### UltramaficContinentalLithosphereStrongZone EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "UltramaficContinentalLithosphereStrongZone"
    )
Markdown.parse(table_str)
```

### UltramaficAsthenosphereDryFertile EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "UltramaficAsthenosphereDryFertile"
    )
Markdown.parse(table_str)
```

### UltramaficAsthenosphereWetFertile EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "UltramaficAsthenosphereWetFertile"
    )
Markdown.parse(table_str)
```

### OceanicGabbroicCrust EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "OceanicGabbroicCrust"
    )
Markdown.parse(table_str)
```

### GabbroicMagma EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "GabbroicMagma"
    )
Markdown.parse(table_str)
```

### LayeredGabbroicCrust EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "LayeredGabbroicCrust"
    )
Markdown.parse(table_str)
```

### LayeredGabbroicMagma EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "LayeredGabbroicMagma"
    )
Markdown.parse(table_str)
```

### SerpentinizedPeridotite EB1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1,
    target_collection_name = "Lithospheric_Deformation_eb1", material_name = "SerpentinizedPeridotite"
    )
Markdown.parse(table_str)
```


