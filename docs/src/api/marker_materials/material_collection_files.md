# Material Collection Files

Material collection files are yamal formatted files that contain groups of materials 
and associated properties. The user can create a custom material collection file or use one of the 
collections that are packaged with EarthBox [MaterialLibrary](@ref).

The general structure of a material collection file is as follows:

```
description: |
  
  Text description goes here.
  - And can include bullets
  - Like these

  And Sections
  ============
  With text that provides additional details....
nmats: 3
MaterialGroupNameA:
    UniqueMaterialNameA:
      property_name_1: [value, units, description]
      property_name_2: [value, units, description]
      ...
    UniqueMaterialNameB
      property_name_1: [value, units, description]
      property__name_2: [value, units, description]
      ...
MaterialGroupNameB
    UniqueMaterialNameC:
      property_name_1: [value, units, description]
      property_name_2: [value, units, description]
      property_name_3: [value, units, description]
      ...
    UniqueMaterialNameD
      property_name_1: [value, units, description]
      property_name_2: [value, units, description]
      ...
```

See the list of valid material property names and standard units for material library files in 
[Material Properties](@ref).

!!! warning "material names"
    All material names in the material collection file must be unique. 
  
!!! info "group names"
    Group names are only for organizational purposes and are ignored by the code.

!!! warning "material property units"
    The units provided in the material library file must match the master units as described in
    [Material Properties](@ref). Material property unit conversion will be added in a future
    version of EarthBox.