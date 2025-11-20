module ArrayCollection

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import EarthBox.Arrays.ArrayTypes.TopoArray2D: TopoArray2DState

"""
    Arrays <: AbstractArrayCollection

Data structure containing array groups for topography evolution.

# Fields
- `gridt::`[`TopoArray2DState`](@ref): `(7, toponum)` Topography grid array
    - ROW 1: x-coordinates of topography grid nodes.
    - ROW 2: y-coordinate of topography grid nodes (i.e. elevation).
    - ROW 3: Updated pre-antidiffusion y-coordinate (i.e. elevation) at topography grid nodes.
    - ROW 4: x-component of velocity at topography grid nodes.
    - ROW 5: y-component of velocity at topography grid nodes.
    - ROW 6: Antidiffusion correction.
    - ROW 7: Thickness of extrusive material at each node.

# Nested Dot Access
- `model.topography.arrays.gridt.array`

# Constructor
    Arrays(toponum::Int)::Arrays

Create a new Arrays collection with the given topography grid resolution.

# Arguments
- `toponum::Int`: Number of topography grid points

"""
mutable struct Arrays <: AbstractArrayCollection
    gridt::TopoArray2DState
end

function Arrays(toponum::Int)::Arrays
    units = [
        "1: m",
        "2: m",
        "3: m",
        "4: m/s",
        "5: m/s",
        "6: m",
        "7: m"
    ]
    descriptions = [
        "1: x-coordinates of topography grid nodes.",
        "2: y-coordinate of topography grid nodes (i.e. elevation).",
        "3: Updated pre-antidiffusion y-coordinate (i.e. elevation) at topography grid nodes.",
        "4: x-component of velocity at topography grid nodes.",
        "5: y-component of velocity at topography grid nodes.",
        "6: Antidiffusion correction.",
        "7: Thickness of extrusive material at each node."
    ]
    gridt = TopoArray2DState(
        initialize_topo_array(toponum),
        "gridt",
        units,
        descriptions
    )
    return Arrays(gridt)
end

function initialize_topo_array(toponum::Int64)::Matrix{Float64}
    return zeros(Float64, 7, toponum)
end

end # module 