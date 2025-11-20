# MaterialLibrary

The `MaterialLibrary` struct provides access to material collections files that
are packaged with EarthBox. See [Material Collection Files](@ref) for more
information about the format of material collection library files.

Material Collection Tables
- [Material Collection: `lithospheric_deformation_eb1`](@ref)

```@docs
EarthBox.MaterialLibraryCollection.MaterialLibrary
```

```@docs
EarthBox.MaterialLibraryCollection.MaterialCollection
```

```@docs
EarthBox.MaterialLibraryCollection.LithosphericDeformationCollections
```

# Example: Accessing Material Collection Files and Material Names

```julia
using EarthBox

lib = MaterialLibraryCollection.MaterialLibrary();

lith_collections = lib.lithospheric_deformation;

collection = lith_collections.lithospheric_deformation_eb1;
path = collection.path
materials = collection.materials

collection = lith_collections.lithospheric_deformation_naliboff17;
path = collection.path
materials = collection.materials

collection = lith_collections.lithospheric_deformation_naliboff17corrected;
lib_path = collection.path
mat_names = collection.materials

collection = lith_collections.lithosphere_deformation_brune14;
lib_path = collection.path
mat_names = collection.materials

collection = lith_collections.lithosphere_deformation_brune14corrected
lib_path = collection.path
mat_names = collection.materials

path = lib.benchmarks.path
materials = lib.benchmarks.materials

path = lib.gerya2019.path
materials = lib.gerya2019.materials

path = lib.viscoelastoplastic.path
materials = lib.viscoelastoplastic.materials

path = lib.simple_mantle_melting.path
materials = lib.simple_mantle_melting.materials

path = lib.seafloor_spreading.path
materials = lib.seafloor_spreading.materials

```
