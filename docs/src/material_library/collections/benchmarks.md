# Material Collection: `benchmarks`

# Material List: Single-material Benchmarks
- [Viscoelastic_Material](@ref)
- [MaterialForSolidBodyRotation](@ref)
- [Layer1RayleighTaylor](@ref)
- [Layer2RayleighTaylor](@ref)
- [IsoviscousFluidBlankenbach89](@ref)
- [IsoviscousFluid1AroundElasticSlab](@ref)
- [IsoviscousFluid2AroundElasticSlab](@ref)
- [ViscoelasticSlabMaterial1](@ref)
- [ViscoelasticSlabMaterial2](@ref)
- [NonNewtonianMaterialChannelFlowBenchmark](@ref)
- [ConstantViscosityChannel](@ref)
- [ConstantViscosityChannelVariableThermalConductivity](@ref)
- [TempDependentViscosityCouetteFlowBenchmark](@ref)

# Material List: Kaus10 Plasticity Benchmark
- [StickyAirKaus10](@ref)
- [CrustKaus10](@ref)
- [WeakSeedKaus10](@ref)

# Material List: Flexure with Triangular Hole Benchmark
- [FlexTestStickyWater](@ref)
- [FlexTestIsoviscousCrust](@ref)
- [FlexTestIsoviscousMantle](@ref)

# Material List: Visco-elastic Extension and Contraction Benchmarks
- [ViscoelasticExtensionStickyAir](@ref)
- [ViscoelasticExtensionStickyWater](@ref)
- [ViscoelasticExtensionSediment](@ref)
- [ViscoelasticExtensionUpperCrust](@ref)
- [ViscoelasticExtensionLowerCrust](@ref)
- [ViscoelasticExtensionUpperCrustStrong](@ref)
- [ViscoelasticExtensionLowerCrustStrong](@ref)
- [ViscoelasticExtensionMantle](@ref)

# Material List: Simple Sedimentation Benchmark
- [SimpleSedStickyAir](@ref)
- [SimpleSedStickyWater](@ref)
- [SimpleSedSedimentWetQuartziteGT95](@ref)
- [SimpleSedPeridotiteDryKK08](@ref)
- [SimpleSedSalt](@ref)

# Material Summary

**Collection File Path**
```
src/material_library/registries/benchmarks/benchmarks.yml
```

# Material Tables

## Material Tables: `Viscoelastic_Stress_Buildup`

### Viscoelastic_Material
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Stress_Buildup", 
    material_name = "Viscoelastic_Material"
    )
Markdown.parse(table_str)
```

## Material Tables: `Solid_Body_Rotation_Benchmark`

### MaterialForSolidBodyRotation
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Solid_Body_Rotation_Benchmark", 
    material_name = "MaterialForSolidBodyRotation"
    )
Markdown.parse(table_str)
```

## Material Tables: `Rayleigh_Taylor_Instability_Benchmark`

### Layer1RayleighTaylor
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Rayleigh_Taylor_Instability_Benchmark", 
    material_name = "Layer1RayleighTaylor"
    )
Markdown.parse(table_str)
```

### Layer2RayleighTaylor
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Rayleigh_Taylor_Instability_Benchmark", 
    material_name = "Layer2RayleighTaylor"
    )
Markdown.parse(table_str)
```

## Material Tables: `Box_Convection_Benchmark`

### IsoviscousFluidBlankenbach89
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Box_Convection_Benchmark", 
    material_name = "IsoviscousFluidBlankenbach89"
    )
Markdown.parse(table_str)
```

## Material Tables: `Elastic_Slab_Benchmark`

### IsoviscousFluid1AroundElasticSlab
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Elastic_Slab_Benchmark", 
    material_name = "IsoviscousFluid1AroundElasticSlab"
    )
Markdown.parse(table_str)
```

### IsoviscousFluid2AroundElasticSlab
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Elastic_Slab_Benchmark", 
    material_name = "IsoviscousFluid2AroundElasticSlab"
    )
Markdown.parse(table_str)
```

### ViscoelasticSlabMaterial1
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Elastic_Slab_Benchmark", 
    material_name = "ViscoelasticSlabMaterial1"
    )
Markdown.parse(table_str)
```

### ViscoelasticSlabMaterial2
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Elastic_Slab_Benchmark", 
    material_name = "ViscoelasticSlabMaterial2"
    )
Markdown.parse(table_str)
```

## Material Tables: `Plasticity_Benchmark_Kaus10`

### StickyAirKaus10
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Plasticity_Benchmark_Kaus10", 
    material_name = "StickyAirKaus10"
    )
Markdown.parse(table_str)
```

### CrustKaus10
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Plasticity_Benchmark_Kaus10", material_name = "CrustKaus10"
    )
Markdown.parse(table_str)
```

### WeakSeedKaus10
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Plasticity_Benchmark_Kaus10", 
    material_name = "WeakSeedKaus10"
    )
Markdown.parse(table_str)
```

## Material Tables: `NonNewtonian_Channel_Flow_Benchmark`

### NonNewtonianMaterialChannelFlowBenchmark
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "NonNewtonian_Channel_Flow_Benchmark", 
    material_name = "NonNewtonianMaterialChannelFlowBenchmark"
    )
Markdown.parse(table_str)
```

## Material Table: `Channel_Flow_Non_Steady_Temperature_Benchmark`

### ConstantViscosityChannel
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Channel_Flow_Non_Steady_Temperature_Benchmark", 
    material_name = "ConstantViscosityChannel"
    )
Markdown.parse(table_str)
```

## Material Tables: `Channel_Flow_Variable_Conductivity`

### ConstantViscosityChannelVariableThermalConductivity
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Channel_Flow_Variable_Conductivity", 
    material_name = "ConstantViscosityChannelVariableThermalConductivity"
    )
Markdown.parse(table_str)
```

# Material Tables: `TempDependentViscosityCouetteFlowBenchmark`

### TempDependentViscosityCouetteFlowBenchmark
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "TempDependentViscosityCouetteFlowBenchmark", 
    material_name = "TempDependentViscosityCouetteFlowBenchmark"
    )
Markdown.parse(table_str)
```

## Material Tables: `FlexureTest_Triangular_Load`

### FlexTestStickyWater
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "FlexureTest_Triangular_Load", 
    material_name = "FlexTestStickyWater"
    )
Markdown.parse(table_str)
```

### FlexTestIsoviscousCrust
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "FlexureTest_Triangular_Load", 
    material_name = "FlexTestIsoviscousCrust"
    )
Markdown.parse(table_str)
```

### FlexTestIsoviscousMantle
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "FlexureTest_Triangular_Load", 
    material_name = "FlexTestIsoviscousMantle"
    )
Markdown.parse(table_str)
```
`
## Material Table: `Viscoelastic_Extension`

### ViscoelasticExtensionStickyAir
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Extension", 
    material_name = "ViscoelasticExtensionStickyAir"
    )
Markdown.parse(table_str)
```

### ViscoelasticExtensionStickyWater
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Extension", 
    material_name = "ViscoelasticExtensionStickyWater"
    )
Markdown.parse(table_str)
```

### ViscoelasticExtensionSediment
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Extension", 
    material_name = "ViscoelasticExtensionSediment"
    )
Markdown.parse(table_str)
```

### ViscoelasticExtensionUpperCrust
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Extension", 
    material_name = "ViscoelasticExtensionUpperCrust"
    )
Markdown.parse(table_str)
```

### ViscoelasticExtensionLowerCrust
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Extension", 
    material_name = "ViscoelasticExtensionLowerCrust"
    )
Markdown.parse(table_str)
```

### ViscoelasticExtensionUpperCrustStrong
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Extension", 
    material_name = "ViscoelasticExtensionUpperCrustStrong"
    )
Markdown.parse(table_str)
```

### ViscoelasticExtensionLowerCrustStrong
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Extension", 
    material_name = "ViscoelasticExtensionLowerCrustStrong"
    )
Markdown.parse(table_str)
```

### ViscoelasticExtensionMantle
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Viscoelastic_Extension", 
    material_name = "ViscoelasticExtensionMantle"
    )
Markdown.parse(table_str)
```

## Material Table: `Simple_Sedimentation`

### SimpleSedStickyAir
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Simple_Sedimentation", 
    material_name = "SimpleSedStickyAir"
    )
Markdown.parse(table_str)
```

### SimpleSedStickyWater
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Simple_Sedimentation", 
    material_name = "SimpleSedStickyWater"
    )
Markdown.parse(table_str)
```

### SimpleSedSedimentWetQuartziteGT95
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Simple_Sedimentation", 
    material_name = "SimpleSedSedimentWetQuartziteGT95"
    )
Markdown.parse(table_str)
```

### SimpleSedPeridotiteDryKK08
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Simple_Sedimentation", 
    material_name = "SimpleSedPeridotiteDryKK08"
    )
Markdown.parse(table_str)
```

### SimpleSedSalt
```@eval
using EarthBox
import Markdown
table_str = make_material_table(
    MaterialLibrary().benchmarks,
    target_collection_name = "Simple_Sedimentation", 
    material_name = "SimpleSedSalt"
    )
Markdown.parse(table_str)
```