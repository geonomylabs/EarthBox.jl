module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .StickyAirGeometryGroup: StickyAirGeometry
import .EarthLayeringGroup: EarthLayering
import .WeakSeedGroup: WeakSeed
import .WeakFaultGroup: WeakFault
import .LithoStrongZonesGroup: LithoStrongZones
import .PlumeGroup: Plume
import .CrustalHoleGroup: CrustalHole
import .RayleighTaylorGroup: RayleighTaylor
import .MobileWallGroup: MobileWall
import .SandboxGroup: Sandbox
import .InternalVelocityZoneGroup: InternalVelocityZone
import .FractureZoneGroup: FractureZone

"""
    Parameters <: AbstractParameterCollection

Collection of geometry parameters.

# Fields
- `sticky_air_geometry::`[`StickyAirGeometry`](@ref): Sticky air layer geometry parameters
- `earth_layering::`[`EarthLayering`](@ref): Earth layering parameters
- `weak_seed::`[`WeakSeed`](@ref): Weak seed parameters
- `weak_fault::`[`WeakFault`](@ref): Weak fault parameters  
- `litho_strong_zones::`[`LithoStrongZones`](@ref): Lithospheric strong zone parameters
- `plume::`[`Plume`](@ref): Mantle plume parameters
- `crustal_hole::`[`CrustalHole`](@ref): Crustal hole parameters
- `rayleigh_taylor::`[`RayleighTaylor`](@ref): Rayleigh-Taylor instability parameters
- `mobile_wall::`[`MobileWall`](@ref): Mobile wall parameters
- `sandbox::`[`Sandbox`](@ref): Sandbox experiment parameters
- `internal_velocity_zone::`[`InternalVelocityZone`](@ref): Internal velocity zone parameters
- `fracture_zone::`[`FractureZone`](@ref): Fracture zone parameters
"""
mutable struct Parameters <: AbstractParameterCollection
    sticky_air_geometry::StickyAirGeometry
    earth_layering::EarthLayering
    weak_seed::WeakSeed
    weak_fault::WeakFault
    litho_strong_zones::LithoStrongZones
    plume::Plume
    crustal_hole::CrustalHole
    rayleigh_taylor::RayleighTaylor
    mobile_wall::MobileWall
    sandbox::Sandbox
    internal_velocity_zone::InternalVelocityZone
    fracture_zone::FractureZone
end

function Parameters()::Parameters
    return Parameters(
        StickyAirGeometry(), # sticky_air_geometry
        EarthLayering(), # earth_layering
        WeakSeed(), # weak_seed
        WeakFault(), # weak_fault
        LithoStrongZones(), # litho_strong_zones
        Plume(), # plume
        CrustalHole(), # crustal_hole
        RayleighTaylor(), # rayleigh_taylor
        MobileWall(), # mobile_wall
        Sandbox(), # sandbox
        InternalVelocityZone(), # internal_velocity_zone
        FractureZone() # fracture_zone
    )
end

end # module 