"""
Module for internal boundary condition arrays.

Provides data structures for specifying internal boundary conditions on velocity
nodes within the computational domain.
"""
module BcInternalGroup

import EarthBox.Arrays.ArrayTypes.InternalBcArrayInt: InternalBcArrayIntState
import EarthBox.Arrays.ArrayTypes.InternalBcArrayFloat: InternalBcArrayFloatState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.bcs.arrays"
const GRP_NAME = "internal"

"""
    BcInternal <: AbstractArrayGroup

Array group for internal boundary conditions.

# Fields
- `bintern_zone::`[`InternalBcArrayIntState`](@ref): (8) : Integer array with grid index information 
    defining spatial zones where internal boundary conditions are applied.
    - (1) Horizontal position of vx nodes with prescribed velocity [no such condition if negative]
    - (2) Min vertical position
    - (3) Max vertical position
    - (4) (not used)
    - (5) Horizontal position of vy nodes with prescribed velocity [no such condition if negative]
    - (6) Min vertical position
    - (7) Max vertical position
    - (8) (not used)
- `bintern_velocity::`[`InternalBcArrayFloatState`](@ref): (8) : Float array containing
    prescribed velocity values for internal boundary conditions.
    - (1) (not used)
    - (2) (not used)
    - (3) (not used)
    - (4) Prescribed shortening velocity, m/s
    - (5) (not used)
    - (6) (not used)
    - (7) (not used)
    - (8) Prescribed vertical velocity, m/s

# Nested Dot Access
- `bintern_zone = $(ROOT_NAME).$(GRP_NAME).bintern_zone.array`
- `bintern_velocity = $(ROOT_NAME).$(GRP_NAME).bintern_velocity.array`

# Constructor
    BcInternal()

"""
mutable struct BcInternal <: AbstractArrayGroup
    bintern_zone::InternalBcArrayIntState
    bintern_velocity::InternalBcArrayFloatState
end

function BcInternal()::BcInternal
    return BcInternal(
        InternalBcArrayIntState(
            fill(-1, 8), # array
            "bintern_zone", # name
            [
                "1: Horizontal position of vx nodes with prescribed velocity [no such condition if negative]",
                "2: Min vertical position",
                "3: Max vertical position",
                "4: (not used)",
                "5: Horizontal position of vy nodes with prescribed velocity [no such condition if negative]",
                "6: Min vertical position",
                "7: Max vertical position",
                "8: (not used)"
            ]
        ),
        InternalBcArrayFloatState(
            zeros(8), # array
            "bintern_velocity", # name
            [
                "1: (not used)",
                "2: (not used)",
                "3: (not used)",
                "4: Prescribed x-velocity, m/s",
                "5: (not used)",
                "6: (not used)",
                "7: (not used)",
                "8: Prescribed y-velocity, m/s"
            ]
        )
    )
end

end # module 