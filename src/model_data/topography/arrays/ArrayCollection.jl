module ArrayCollection

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import EarthBox.Arrays.ArrayTypes.TopoArray2D: TopoArray2DState
import EarthBox.Arrays.ArrayTypes.Array3DFloat: Array3DFloatState
import EarthBox.Arrays.ArrayTypes.Array1DFloat: Array1DFloatState

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
- `layer_tops_buffer::`[`Array1DFloatState`](@ref): `(toponum)` Pre-allocated
    scratch buffer used by `TopAndBottom.calculate_top_and_bottom_of_layer_opt`
    via the `tops_buffer` keyword to avoid allocating a fresh `tops` vector
    each call. Caller zeros it before use.
- `layer_bottoms_buffer::`[`Array1DFloatState`](@ref): `(toponum)` Same role
    as `layer_tops_buffer` but for the `bottoms` array.
- `oceanic_moho_buffer::`[`Array1DFloatState`](@ref): `(toponum)` Persistent
    output buffer for the smoothed oceanic moho y-coordinates produced by
    `Fractionation.calculate_oceanic_moho`. Reused each fractionation call.
- `partial_melt_buffer::`[`Array1DFloatState`](@ref): `(toponum)` Persistent
    output buffer for the smoothed top-of-mantle-partial-melt y-coordinates
    produced by `Drainage.calculate_top_of_mantle_partial_melt_domain`.
    Reused each fractionation/drainage call.

# Nested Dot Access
- `model.topography.arrays.gridt.array`
- `model.topography.arrays.compaction_array.array`
- `model.topography.arrays.layer_tops_buffer.array`
- `model.topography.arrays.layer_bottoms_buffer.array`
- `model.topography.arrays.oceanic_moho_buffer.array`
- `model.topography.arrays.partial_melt_buffer.array`

# Constructor
    Arrays(toponum::Int)::Arrays

Create a new Arrays collection with the given topography grid resolution.

# Arguments
- `toponum::Int`: Number of topography grid points

"""
mutable struct Arrays <: AbstractArrayCollection
    gridt::TopoArray2DState
    compaction_array::Array3DFloatState
    layer_tops_buffer::Array1DFloatState
    layer_bottoms_buffer::Array1DFloatState
    oceanic_moho_buffer::Array1DFloatState
    partial_melt_buffer::Array1DFloatState
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
    layer_tops_buffer = Array1DFloatState(
        zeros(Float64, toponum),
        "layer_tops_buffer",
        "m",
        "`(toponum)` : Pre-allocated scratch buffer for the `tops` output of "
        * "TopAndBottom.calculate_top_and_bottom_of_layer_opt. Reused across "
        * "fractionation and drainage calls."
    )
    layer_bottoms_buffer = Array1DFloatState(
        zeros(Float64, toponum),
        "layer_bottoms_buffer",
        "m",
        "`(toponum)` : Pre-allocated scratch buffer for the `bottoms` output "
        * "of TopAndBottom.calculate_top_and_bottom_of_layer_opt. Reused "
        * "across fractionation and drainage calls."
    )
    oceanic_moho_buffer = Array1DFloatState(
        zeros(Float64, toponum),
        "oceanic_moho_buffer",
        "m",
        "`(toponum)` : Persistent output buffer for smoothed oceanic moho "
        * "y-coordinates produced by Fractionation.calculate_oceanic_moho."
    )
    partial_melt_buffer = Array1DFloatState(
        zeros(Float64, toponum),
        "partial_melt_buffer",
        "m",
        "`(toponum)` : Persistent output buffer for smoothed top-of-"
        * "mantle-partial-melt y-coordinates produced by "
        * "Drainage.calculate_top_of_mantle_partial_melt_domain."
    )
    return Arrays(
        gridt, compaction_array,
        layer_tops_buffer, layer_bottoms_buffer,
        oceanic_moho_buffer, partial_melt_buffer
    )
end

function initialize_topo_array(toponum::Int64)::Matrix{Float64}
    return zeros(Float64, 7, toponum)
end

end # module
