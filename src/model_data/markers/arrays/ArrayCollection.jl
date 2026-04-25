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
import .AdvectionGroup: Advection
import .StructureGroup: Structure
import .CompactionGroup: Compaction
import .SolidificationGroup: Solidification
import .RecycleGroup: Recycle
import .LithostaticScratchGroup: LithostaticScratch
import .SerpentinizationGroup: Serpentinization
import .SubgridHeatGroup: SubgridHeat

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
- `advection::`[`Advection`](@ref): Marker advection velocity/spin arrays
- `structure::`[`Structure`](@ref): Marker structure-finding scratch buffers
- `compaction::`[`Compaction`](@ref): Marker compaction scratch buffers
- `solidification::`[`Solidification`](@ref): Marker solidification scratch buffers
- `recycle::`[`Recycle`](@ref): Marker recycling scratch buffers
- `lithostatic_scratch::`[`LithostaticScratch`](@ref): Lithostatic-pressure column-filter scratch
- `serpentinization::`[`Serpentinization`](@ref): Serpentinization scratch buffers
- `subgrid_heat::`[`SubgridHeat`](@ref): Subgrid heat-diffusion scratch buffers
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
    advection::Advection
    structure::Structure
    compaction::Compaction
    solidification::Solidification
    recycle::Recycle
    lithostatic_scratch::LithostaticScratch
    serpentinization::Serpentinization
    subgrid_heat::SubgridHeat
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
        Melt(marknum),
        Advection(marknum),
        Structure(marknum),
        Compaction(marknum),
        Solidification(marknum),
        Recycle(marknum),
        LithostaticScratch(marknum),
        Serpentinization(marknum),
        SubgridHeat(marknum)
    )
end

end # module
