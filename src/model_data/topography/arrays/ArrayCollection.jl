module ArrayCollection

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import EarthBox.Arrays.ArrayTypes.TopoArray2D: TopoArray2DState
import EarthBox.Arrays.ArrayTypes.Array3DFloat: Array3DFloatState

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
- `compaction_array::`[`Array3DFloatState`](@ref): `(toponum, 20, 9)` Pre-allocated
    scratch buffer used by `MarkerCompaction.compact_sediment_and_advect_markers!`.
    Caller zeros it via `fill!` at start of each call. Only valid during a
    single compaction invocation; not persistent state.

# Nested Dot Access
- `model.topography.arrays.gridt.array`
- `model.topography.arrays.compaction_array.array`

# Constructor
    Arrays(toponum::Int)::Arrays

Create a new Arrays collection with the given topography grid resolution.

# Arguments
- `toponum::Int`: Number of topography grid points

"""
mutable struct Arrays <: AbstractArrayCollection
    gridt::TopoArray2DState
    compaction_array::Array3DFloatState
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
    compaction_array = Array3DFloatState(
        zeros(Float64, toponum, 20, 9),
        "compaction_array",
        "mixed",
        "`(toponum, 20, 9)` : Pre-allocated scratch buffer for " *
        "MarkerCompaction.compact_sediment_and_advect_markers!. The caller " *
        "zeros it via fill! at start of each call; contents are not persistent " *
        "state. Third-dim slots: 1=y-coord, 2=initial porosity, 3=decay depth, " *
        "4=marker count, 5=thickness, 6=thickness delta, 7=cumulative " *
        "y-displacement, 8=max burial depth, 9=updated burial depth."
    )
    return Arrays(gridt, compaction_array)
end

function initialize_topo_array(toponum::Int64)::Matrix{Float64}
    return zeros(Float64, 7, toponum)
end

end # module
