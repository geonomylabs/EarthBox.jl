"""
Module for subgrid heat diffusion scratch buffers.

Single marker-sized buffer reused by
`SubgridDiffusion.calculate_subgrid_temperature_change_and_correct_marker_temperature!`
to avoid a marknum-sized `Vector{Float64}` allocation per heat-solver call.
"""
module SubgridHeatGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "subgrid_heat"

const ADATA = get_eb_arrays()

"""
    SubgridHeat <: AbstractArrayGroup

Array group holding scratch buffers reused by subgrid heat-diffusion routines.

# Fields
- `marker_subgrid_temp_delta_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64

# Constructor
    SubgridHeat(marknum::Int)
"""
mutable struct SubgridHeat <: AbstractArrayGroup
    marker_subgrid_temp_delta_buffer::MarkerArrayFloat1DState{Float64}
end

function SubgridHeat(marknum::Int)::SubgridHeat
    return SubgridHeat(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_subgrid_temp_delta_buffer.name,
            ADATA.marker_subgrid_temp_delta_buffer.units,
            ADATA.marker_subgrid_temp_delta_buffer.description
        )
    )
end

end # module
