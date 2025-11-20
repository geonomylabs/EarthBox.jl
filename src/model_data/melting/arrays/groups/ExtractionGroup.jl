module ExtractionGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.melting.arrays"
const GRP_NAME = "extraction"

const ADATA = get_eb_arrays()

"""
    Extraction <: AbstractArrayGroup

Array group containing drainage basin tracking arrays for melt extraction.

# Fields
- `xstart_drainage::`[`GridArray1DState`](@ref): $(ADATA.xstart_drainage.description)
- `xend_drainage::`[`GridArray1DState`](@ref): $(ADATA.xend_drainage.description)
- `melt_residuals::`[`GridArray1DState`](@ref): $(ADATA.melt_residuals.description)
- `extrusion_volumes::`[`GridArray1DState`](@ref): $(ADATA.extrusion_volumes.description)
- `xmid_molten_zones::`[`GridArray1DState`](@ref): $(ADATA.xmid_molten_zones.description)
- `ytop_molten_zones::`[`GridArray1DState`](@ref): $(ADATA.ytop_molten_zones.description)
- `width_molten_zones::`[`GridArray1DState`](@ref): $(ADATA.width_molten_zones.description)
- `avg_shallow_partial_melt_xcoors::`[`GridArray1DState`](@ref): $(ADATA.avg_shallow_partial_melt_xcoors.description)
- `avg_shallow_partial_melt_ycoors::`[`GridArray1DState`](@ref): $(ADATA.avg_shallow_partial_melt_ycoors.description)
- `xstart_drainage_o::`[`GridArray1DState`](@ref): $(ADATA.xstart_drainage_o.description)
- `xend_drainage_o::`[`GridArray1DState`](@ref): $(ADATA.xend_drainage_o.description)

# Nested Dot Access
- `xstart_drainage = $(ROOT_NAME).$(GRP_NAME).xstart_drainage.array`
- `xend_drainage = $(ROOT_NAME).$(GRP_NAME).xend_drainage.array`
- `melt_residuals = $(ROOT_NAME).$(GRP_NAME).melt_residuals.array`
- `extrusion_volumes = $(ROOT_NAME).$(GRP_NAME).extrusion_volumes.array`
- `xmid_molten_zones = $(ROOT_NAME).$(GRP_NAME).xmid_molten_zones.array`
- `ytop_molten_zones = $(ROOT_NAME).$(GRP_NAME).ytop_molten_zones.array`
- `width_molten_zones = $(ROOT_NAME).$(GRP_NAME).width_molten_zones.array`
- `avg_shallow_partial_melt_xcoors = $(ROOT_NAME).$(GRP_NAME).avg_shallow_partial_melt_xcoors.array`
- `avg_shallow_partial_melt_ycoors = $(ROOT_NAME).$(GRP_NAME).avg_shallow_partial_melt_ycoors.array`
- `xstart_drainage_o = $(ROOT_NAME).$(GRP_NAME).xstart_drainage_o.array`
- `xend_drainage_o = $(ROOT_NAME).$(GRP_NAME).xend_drainage_o.array`

# Constructor
    Extraction()

Create a new Extraction array group with default drainage basin arrays.

# Returns
- `Extraction`: New Extraction array group with initialized arrays

"""
mutable struct Extraction <: AbstractArrayGroup
    xstart_drainage::GridArray1DState
    xend_drainage::GridArray1DState
    melt_residuals::GridArray1DState
    extrusion_volumes::GridArray1DState
    xmid_molten_zones::GridArray1DState
    ytop_molten_zones::GridArray1DState
    width_molten_zones::GridArray1DState
    avg_shallow_partial_melt_xcoors::GridArray1DState
    avg_shallow_partial_melt_ycoors::GridArray1DState
    xstart_drainage_o::GridArray1DState
    xend_drainage_o::GridArray1DState
end

function Extraction()
    nmax_basin = 100  # maximum number of permitted drainage basins
    return Extraction(
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.xstart_drainage.name,            # name
            ADATA.xstart_drainage.units,           # units
            ADATA.xstart_drainage.description      # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.xend_drainage.name,              # name
            ADATA.xend_drainage.units,             # units
            ADATA.xend_drainage.description        # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.melt_residuals.name,             # name
            ADATA.melt_residuals.units,            # units
            ADATA.melt_residuals.description       # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.extrusion_volumes.name,          # name
            ADATA.extrusion_volumes.units,         # units
            ADATA.extrusion_volumes.description    # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.xmid_molten_zones.name,          # name
            ADATA.xmid_molten_zones.units,         # units
            ADATA.xmid_molten_zones.description    # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.ytop_molten_zones.name,          # name
            ADATA.ytop_molten_zones.units,         # units
            ADATA.ytop_molten_zones.description    # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.width_molten_zones.name,         # name
            ADATA.width_molten_zones.units,        # units
            ADATA.width_molten_zones.description   # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),                    # array
            ADATA.avg_shallow_partial_melt_xcoors.name,    # name
            ADATA.avg_shallow_partial_melt_xcoors.units,   # units
            ADATA.avg_shallow_partial_melt_xcoors.description # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),                    # array
            ADATA.avg_shallow_partial_melt_ycoors.name,    # name
            ADATA.avg_shallow_partial_melt_ycoors.units,   # units
            ADATA.avg_shallow_partial_melt_ycoors.description # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.xstart_drainage_o.name,          # name
            ADATA.xstart_drainage_o.units,         # units
            ADATA.xstart_drainage_o.description    # description
        ),
        GridArray1DState(
            zeros(Float64, nmax_basin),           # array
            ADATA.xend_drainage_o.name,            # name
            ADATA.xend_drainage_o.units,           # units
            ADATA.xend_drainage_o.description      # description
        )
    )
end

end # module
