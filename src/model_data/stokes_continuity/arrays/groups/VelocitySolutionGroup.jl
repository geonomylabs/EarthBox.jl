module VelocitySolutionGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.SolutionArray1D: SolutionArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "velocity_solution"

const ADATA = get_eb_arrays()

"""
    VelocitySolution <: AbstractArrayGroup

Array group containing velocity solution arrays.

# Fields
- `vxy_old::`[`SolutionArray1DState`](@ref): $(ADATA.vxy_old.description)
- `vxy::`[`SolutionArray1DState`](@ref): $(ADATA.vxy.description)

# Nested Dot Access
- `vxy_old = $(ROOT_NAME).$(GRP_NAME).vxy_old.array`
- `vxy = $(ROOT_NAME).$(GRP_NAME).vxy.array`

# Constructor
    VelocitySolution(ynum::Int, xnum::Int)

Create a new VelocitySolution array group with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `VelocitySolution`: New VelocitySolution array group with initialized arrays

"""
mutable struct VelocitySolution <: AbstractArrayGroup
    vxy_old::SolutionArray1DState
    vxy::SolutionArray1DState
end

function VelocitySolution(ynum::Int, xnum::Int)::VelocitySolution
    return VelocitySolution(
        SolutionArray1DState(
            ynum,                        # ynum
            xnum,                        # xnum
            ADATA.vxy_old.name,          # name
            ADATA.vxy_old.units,         # units
            ADATA.vxy_old.description,   # description
            "velocity"                   # array_type
        ),
        SolutionArray1DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.vxy.name,          # name
            ADATA.vxy.units,         # units
            ADATA.vxy.description,   # description
            "velocity"               # array_type
        )
    )
end

end # module
