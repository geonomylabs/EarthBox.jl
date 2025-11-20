module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .GridOptionsGroup: GridOptions
import .GridGeometryGroup: GridGeometry
import .GridRefinementGroup: GridRefinement

"""
    Parameters <: AbstractParameterCollection

Parameter collection for 2D staggered grids.

# Fields
- `grid_options::`[`GridOptions`](@ref): Grid type options
- `geometry::`[`GridGeometry`](@ref): Dimensions and resolution of basic grids
- `refinement::`[`GridRefinement`](@ref): T-type grid refinement parameters

# Constructor
    Parameters(ynum::Int, xnum::Int, ysize::Float64, xsize::Float64)

Create a new Parameters collection with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction  
- `ysize::Float64`: Domain size in y-direction in meters
- `xsize::Float64`: Domain size in x-direction in meters
"""
mutable struct Parameters <: AbstractParameterCollection
    grid_options::GridOptions
    geometry::GridGeometry
    refinement::GridRefinement
end

function Parameters(
    ynum::Int, xnum::Int, ysize::Float64, xsize::Float64
)::Parameters
    return Parameters(
        GridOptions(),
        GridGeometry(ynum, xnum, ysize, xsize),
        GridRefinement()
    )
end

end # module 