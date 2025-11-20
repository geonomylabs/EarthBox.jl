module InitializeInterpolation

include("utils/BilinearDenominator.jl")
include("utils/MarkerGridNodeWeights.jl")

import Profile
using EarthBox.ModelDataContainer: ModelData
using ..MarkerGridMapping
using .MarkerGridNodeWeights
using .BilinearDenominator

"""
    initialize_bilinear_interpolation(model::ModelData)::Nothing

Initialize bilinear interpolation parameters.

This method manages functions that calculate for each marker the indices
and distances of the upper left basic node and weights for all four basic
nodes and central node of the basic cell that contains the marker. These
marker-grid relationships are calculated for the basic grid and the 
vy-grid. The function also calculates the summed marker weights on nodes
of the basic, pressure and vy-grids that are used in the denominator of
the bilinear interpolation equations.

Several arrays that are stored in the main model data container
are updated by functions called by the current method including the
following:

- `marker_yn::Vector{Int32}` (nmarkers): y-index of the upper left basic node for each marker
- `marker_xn::Vector{Int32}` (nmarkers): x-index of the upper left basic node for each marker
- `marker_yn_vy::Vector{Int32}` (nmarkers): y-index of the upper left vy node for each marker
- `marker_xn_vy::Vector{Int32}` (nmarkers): x-index of the upper left vy node for each marker
- `marker_dx::Vector{Float64}` (nmarkers) : normalized x-distance from upper left basic node
- `marker_dy::Vector{Float64}` (nmarkers): normalized y-distance from upper left basic node
- `marker_dx_vy::Vector{Float64}` (nmarkers): normalized x-distance from upper left vy node
- `marker_dy_vy::Vector{Float64}` (nmarkers): normalized y-distance from upper left vy node
- `marker_wtforULnode::Vector{Float64}` (nmarkers): weight for upper left basic node for each marker
- `marker_wtforURnode::Vector{Float64}` (nmarkers): weight for upper right basic node for each marker
- `marker_wtforLLnode::Vector{Float64}` (nmarkers): weight for lower left basic node for each marker
- `marker_wtforLRnode::Vector{Float64}` (nmarkers): weight for lower right basic node for each marker
- `marker_wtforCnode::Vector{Float64}` (nmarkers): weight for center basic node for each marker
- `marker_wtforULnodeVy::Vector{Float64}` (nmarkers): weight for upper left vy node for each marker
- `marker_wtforURnodeVy::Vector{Float64}` (nmarkers): weight for upper right vy node for each marker
- `marker_wtforLLnodeVy::Vector{Float64}` (nmarkers): weight for lower left vy node for each marker
- `marker_wtforLRnodeVy::Vector{Float64}` (nmarkers): weight for lower right vy node for each marker
- `wtnodes::Matrix{Float64}` (ynum, xnum): summed weights for basic nodes using an inclusive approach
- `wtetas::Matrix{Float64}` (ynum, xnum): summed weights for basic nodes using an exclusive approach 
- `wtetan::Matrix{Float64}` (ynum-1, xnum-1): summed weights for pressure nodes using an inclusive approach
- `wtnodes_vy::Matrix{Float64}` (ynum, xnum+1): summed weights for vy nodes using an inclusive approach

Inclusive refers to an approach where all markers within a grid cell are
included in the calculation of the weight sum for each node. Whereas
exclusive refers to an approach where only markers within a search
distance of a node are included in the calculation of the weight sum
for that node.

# Arguments
- `model::ModelData`: Model data container containing arrays and parameters

# Returns
- `nothing`
"""
function initialize_bilinear_interpolation!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    # Marker-grid relationships on basic and pressure grids
    MarkerGridMapping.calc_upper_left_basic_grid_cell_indices_for_markers!(model, inside_flags)
    MarkerGridMapping.calc_normalized_distances_from_upper_left_basic_node!(model, inside_flags)
    MarkerGridNodeWeights.calc_marker_weights_for_basic_and_pressure_grids!(model, inside_flags)
    # Marker-grid relationships on Vy staggered grid
    MarkerGridMapping.calc_upper_left_vy_grid_cell_indices_for_markers!(model, inside_flags)
    MarkerGridMapping.calc_normalized_distances_from_upper_left_vy_node!(model, inside_flags)
    MarkerGridNodeWeights.calc_marker_weights_for_vy_grid!(model, inside_flags)
    # Clear marker weight sums on all grids
    BilinearDenominator.clear_marker_weight_sums_on_grids!(model)
    # Calculate marker weight sums for basic and pressure nodes
    BilinearDenominator.calc_marker_weight_sums_for_basic_and_pressure_nodes!(model, inside_flags)
    # Calculate marker weight sums for Vy nodes
    BilinearDenominator.calc_marker_weight_sums_for_vy_nodes!(model, inside_flags)
    return nothing
end

end # module 