module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .GridOptionsGroup: GridOptions
import .GridGeometryGroup: GridGeometry
import .GridRefinementGroup: GridRefinement

mutable struct Parameters <: AbstractParameterCollection
    grid_options::GridOptions
    geometry::GridGeometry
    refinement::GridRefinement
end

function Parameters(
    ynum::Int64,
    xnum::Int64,
    znum::Int64,
    ysize::Float64,
    xsize::Float64,
    zsize::Float64
)::Parameters
    return Parameters(
        GridOptions(),
        GridGeometry(ynum, xnum, znum, ysize, xsize, zsize),
        GridRefinement()
    )
end

end # module 