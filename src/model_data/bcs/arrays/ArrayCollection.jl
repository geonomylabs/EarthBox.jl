module ArrayCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .TemperatureGroup: Temperature
import .VelocityGroup: Velocity
import .VelCompGroup: VelComp
import .BcInternalGroup: BcInternal

"""
    Arrays

Collection of boundary condition arrays.

# Fields
- `temperature::`[`Temperature`](@ref): Temperature boundary condition arrays
- `velocity::`[`Velocity`](@ref): Velocity boundary condition arrays
- `vel_comp::`[`VelComp`](@ref): Velocity component arrays
- `internal::`[`BcInternal`](@ref): Internal boundary condition arrays
"""
mutable struct Arrays <: AbstractArrayCollection
    temperature::Temperature
    velocity::Velocity
    vel_comp::VelComp
    internal::BcInternal
end

function Arrays(ynum::Int, xnum::Int)::Arrays
    return Arrays(
        Temperature(ynum, xnum),
        Velocity(ynum, xnum),
        VelComp(ynum, xnum),
        BcInternal()
    )
end

end # module 