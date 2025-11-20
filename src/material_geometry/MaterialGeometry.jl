module MaterialGeometry

include("crustal_hole/CrustalHole.jl")
include("earth_layering/EarthLayering.jl")
include("fracture_zone/FractureZone.jl")
include("internal_velocity_zone/InternalVelocityZone.jl")
include("litho_strong_zones/LithoStrongZones.jl")
include("mobile_wall/MobileWall.jl")
include("plume/Plume.jl")
include("rayleigh_taylor/RayleighTaylor.jl")
include("sandbox/Sandbox.jl")
include("sticky_air_geometry/StickyAirGeometry.jl")
include("weak_fault/WeakFault.jl")
include("weak_seed/WeakSeed.jl")

end