module MarkerContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import ..Grids2dContainer: Grids
import .ArrayCollection: Arrays
import .ParameterCollection: Parameters
import ..DocTools: make_collection_doc

export Markers

"""
    Markers <: CollectionContainer

Data structure containing parameter and array collections for markers.

# Fields
- `parameters::`[`Parameters`](@ref): Parameter groups for marker including distribution, advection, 
    diffusion, and recycling settings
- `arrays::`[`Arrays`](@ref): Array groups for marker arrays containing material properties, 
    state variables, and location information

"""
mutable struct Markers <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

"""
    Markers(marker_parameters::NamedTuple)::Markers

Construct a new Markers struct with the given marker parameters.

# Arguments
- `marker_parameters`: A named tuple containing marker distribution parameters.

The named tuple `marker_parameters` must contain the following parameters:
- `dx_marker`: Marker spacing in x-direction in meters
- `dy_marker`: Marker spacing in y-direction in meters
- `nmarkers_cell_x`: Number of markers per cell in x-direction
- `nmarkers_cell_y`: Number of markers per cell in y-direction
- `mynum`: Number of markers in y-direction
- `mxnum`: Number of markers in x-direction
- `marknum`: Total number of markers
- `mystep`: Marker spacing in y-direction
- `mxstep`: Marker spacing in x-direction

# Returns
- A new Markers struct with the given marker parameters
"""
function Markers(marker_parameters::NamedTuple)::Markers
    @assert marker_parameters.nmarkers_cell_x >= 1.0 "Number of markers per cell in x-direction must be greater than or equal to 1"
    @assert marker_parameters.nmarkers_cell_y >= 1.0 "Number of markers per cell in y-direction must be greater than or equal to 1"
    parameters = Parameters(marker_parameters)
    marknum = parameters.distribution.marknum.value
    @assert marknum > 0 "Total number of markers must be greater than 0"
    return Markers(parameters, Arrays(marknum))
end

end # module 