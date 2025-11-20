module PlasticDefGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.stokes_continuity.arrays"
const GRP_NAME = "plastic_def"

const ADATA = get_eb_arrays()

"""
    PlasticDef <: AbstractArrayGroup

Array group containing plastic deformation arrays.

# Fields
- `plastics::`[`ScalarArray2DState`](@ref): $(ADATA.plastics.description)
- `plastics0::`[`ScalarArray2DState`](@ref): $(ADATA.plastics0.description)
- `plasticn::`[`ScalarArray2DState`](@ref): $(ADATA.plasticn.description)
- `plasticn0::`[`ScalarArray2DState`](@ref): $(ADATA.plasticn0.description)
- `plastic_yield::`[`ScalarArray2DState`](@ref): $(ADATA.plastic_yield.description)
- `yield_error::`[`ScalarArray2DState`](@ref): $(ADATA.yield_error.description)
- `cohesion_grid::`[`ScalarArray2DState`](@ref): $(ADATA.cohesion_grid.description)
- `cohesion_grid0::`[`ScalarArray2DState`](@ref): $(ADATA.cohesion_grid0.description)
- `fric_degrees_grid::`[`ScalarArray2DState`](@ref): $(ADATA.fric_degrees_grid.description)
- `fric_degrees_grid0::`[`ScalarArray2DState`](@ref): $(ADATA.fric_degrees_grid0.description)
- `dilatation_grid::`[`ScalarArray2DState`](@ref): $(ADATA.dilatation_grid.description)
- `dilatation_grid0::`[`ScalarArray2DState`](@ref): $(ADATA.dilatation_grid0.description)

# Nested Dot Access
- `plastics = $(ROOT_NAME).$(GRP_NAME).plastics.array`
- `cohesion_grid = $(ROOT_NAME).$(GRP_NAME).cohesion_grid.array`

# Constructor
    PlasticDef(ynum::Int, xnum::Int)

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

# Returns
- `PlasticDef`: New PlasticDef array group with initialized arrays

"""
mutable struct PlasticDef <: AbstractArrayGroup
    plastics::ScalarArray2DState
    plastics0::ScalarArray2DState
    plasticn::ScalarArray2DState
    plasticn0::ScalarArray2DState
    plastic_yield::ScalarArray2DState
    yield_error::ScalarArray2DState
    cohesion_grid::ScalarArray2DState
    cohesion_grid0::ScalarArray2DState
    fric_degrees_grid::ScalarArray2DState
    fric_degrees_grid0::ScalarArray2DState
    dilatation_grid::ScalarArray2DState
    dilatation_grid0::ScalarArray2DState
end

function PlasticDef(ynum::Int, xnum::Int)::PlasticDef
    data = PlasticDef(
        ScalarArray2DState(
            ynum,                        # ynum
            xnum,                        # xnum
            ADATA.plastics.name,         # name
            ADATA.plastics.units,        # units
            ADATA.plastics.grid_type,    # grid_type
            ADATA.plastics.description   # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.plastics0.name,         # name
            ADATA.plastics0.units,        # units
            ADATA.plastics0.grid_type,    # grid_type
            ADATA.plastics0.description   # description
        ),
        ScalarArray2DState(
            ynum,                        # ynum
            xnum,                        # xnum
            ADATA.plasticn.name,         # name
            ADATA.plasticn.units,        # units
            ADATA.plasticn.grid_type,    # grid_type
            ADATA.plasticn.description   # description
        ),
        ScalarArray2DState(
            ynum,                         # ynum
            xnum,                         # xnum
            ADATA.plasticn0.name,         # name
            ADATA.plasticn0.units,        # units
            ADATA.plasticn0.grid_type,    # grid_type
            ADATA.plasticn0.description   # description
        ),
        ScalarArray2DState(
            ynum,                           # ynum
            xnum,                           # xnum
            ADATA.plastic_yield.name,       # name
            ADATA.plastic_yield.units,      # units
            ADATA.plastic_yield.grid_type,  # grid_type
            ADATA.plastic_yield.description # description
        ),
        ScalarArray2DState(
            ynum,                        # ynum
            xnum,                        # xnum
            ADATA.yield_error.name,      # name
            ADATA.yield_error.units,     # units
            ADATA.yield_error.grid_type, # grid_type
            ADATA.yield_error.description # description
        ),
        ScalarArray2DState(
            ynum,                           # ynum
            xnum,                           # xnum
            ADATA.cohesion_grid.name,       # name
            ADATA.cohesion_grid.units,      # units
            ADATA.cohesion_grid.grid_type,  # grid_type
            ADATA.cohesion_grid.description # description
        ),
        ScalarArray2DState(
            ynum,                             # ynum
            xnum,                             # xnum
            ADATA.cohesion_grid0.name,        # name
            ADATA.cohesion_grid0.units,       # units
            ADATA.cohesion_grid0.grid_type,   # grid_type
            ADATA.cohesion_grid0.description  # description
        ),
        ScalarArray2DState(
            ynum,                                # ynum
            xnum,                                # xnum
            ADATA.fric_degrees_grid.name,        # name
            ADATA.fric_degrees_grid.units,       # units
            ADATA.fric_degrees_grid.grid_type,   # grid_type
            ADATA.fric_degrees_grid.description  # description
        ),
        ScalarArray2DState(
            ynum,                                 # ynum
            xnum,                                 # xnum
            ADATA.fric_degrees_grid0.name,        # name
            ADATA.fric_degrees_grid0.units,       # units
            ADATA.fric_degrees_grid0.grid_type,   # grid_type
            ADATA.fric_degrees_grid0.description  # description
        ),
        ScalarArray2DState(
            ynum,                              # ynum
            xnum,                              # xnum
            ADATA.dilatation_grid.name,        # name
            ADATA.dilatation_grid.units,       # units
            ADATA.dilatation_grid.grid_type,   # grid_type
            ADATA.dilatation_grid.description  # description
        ),
        ScalarArray2DState(
            ynum,                               # ynum
            xnum,                               # xnum
            ADATA.dilatation_grid0.name,        # name
            ADATA.dilatation_grid0.units,       # units
            ADATA.dilatation_grid0.grid_type,   # grid_type
            ADATA.dilatation_grid0.description  # description
        )
    )
    update_output_format(data)
    return data
end

function update_output_format(data::PlasticDef)::Nothing
    data.plastics0.outform.fname = "plastics0"
    data.plastics.outform.fname = "plastics"
    data.plasticn0.outform.fname = "plasticn0"
    data.plasticn.outform.fname = "plasticn"
    data.cohesion_grid.outform.fname = "cohesion_grid_Pa"
    data.fric_degrees_grid.outform.fname = "fric_degrees_grid"
    data.cohesion_grid0.outform.fname = "cohesion_grid0_Pa"
    data.fric_degrees_grid0.outform.fname = "fric_degrees_grid0"
    data.dilatation_grid.outform.fname = "dilatation_grid_degrees"
    return nothing
end

end # module
