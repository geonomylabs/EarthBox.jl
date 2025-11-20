module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .StratigraphyGroup: Stratigraphy
import .RheologyGroup: Rheology
import .ThermalGroup: Thermal
import .MaterialGroup: Material
import .StressGroup: Stress
import .StrainGroup: Strain
import .PressureGroup: Pressure
import .LocationGroup: Location
import .GridMarkerRelationshipGroup: GridMarkerRelationship
import .MeltGroup: Melt

"""
    Arrays <: AbstractArrayCollection

Collection of marker arrays.

# Fields
- `strat::`[`Stratigraphy`](@ref): Stratigraphy marker arrays
- `rheology::`[`Rheology`](@ref): Rheology marker arrays
- `thermal::`[`Thermal`](@ref): Thermal marker arrays
- `material::`[`Material`](@ref): Material marker arrays
- `stress::`[`Stress`](@ref): Stress marker arrays
- `strain::`[`Strain`](@ref): Strain marker arrays
- `pressure::`[`Pressure`](@ref): Pressure marker arrays
- `location::`[`Location`](@ref): Location marker arrays
- `grid_marker_relationship::`[`GridMarkerRelationship`](@ref): Grid-marker relationship arrays
- `melt::`[`Melt`](@ref): Melt marker arrays
"""
mutable struct Arrays <: AbstractArrayCollection
    strat::Stratigraphy
    rheology::Rheology
    thermal::Thermal
    material::Material
    stress::Stress
    strain::Strain
    pressure::Pressure
    location::Location
    grid_marker_relationship::GridMarkerRelationship
    melt::Melt
end

function Arrays(marknum::Int)::Arrays
    return Arrays(
        Stratigraphy(marknum),
        Rheology(marknum),
        Thermal(marknum),
        Material(marknum),
        Stress(marknum),
        Strain(marknum),
        Pressure(marknum),
        Location(marknum),
        GridMarkerRelationship(marknum),
        Melt(marknum)
    )
end

end # module
