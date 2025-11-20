module StokesSolutionGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.SolutionArray1D: SolutionArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "stokes_solution"

const ADATA = get_eb_arrays()

"""
    StokesSolution <: AbstractArrayGroup

Array group containing combined Stokes solution arrays.

# Fields
- `soluv1_old::`[`SolutionArray1DState`](@ref): $(ADATA.soluv1_old.description)
- `soluv1::`[`SolutionArray1DState`](@ref): $(ADATA.soluv1.description)

# Nested Dot Access
- `soluv1_old = $(ROOT_NAME).$(GRP_NAME).soluv1_old.array`
- `soluv1 = $(ROOT_NAME).$(GRP_NAME).soluv1.array`

# Constructor
    StokesSolution(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `StokesSolution`: New StokesSolution array group with initialized arrays

"""
mutable struct StokesSolution <: AbstractArrayGroup
    soluv1_old::SolutionArray1DState
    soluv1::SolutionArray1DState
end

function StokesSolution(ynum::Int, xnum::Int)::StokesSolution
    return StokesSolution(
        SolutionArray1DState(
            ynum,                        # ynum
            xnum,                        # xnum
            ADATA.soluv1_old.name,       # name
            ADATA.soluv1_old.units,      # units
            ADATA.soluv1_old.description, # description
            "normal"                     # array_type
        ),
        SolutionArray1DState(
            ynum,                    # ynum
            xnum,                    # xnum
            ADATA.soluv1.name,       # name
            ADATA.soluv1.units,      # units
            ADATA.soluv1.description, # description
            "normal"                 # array_type
        )
    )
end

end # module
