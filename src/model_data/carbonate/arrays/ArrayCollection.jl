module ArrayCollection

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import EarthBox.Arrays.ArrayTypes.CarbArray2D: CarbArray2DState

const ROOT_NAME = "model.carbonate.arrays"

#const ADATA = get_eb_arrays()

"""
    Arrays <: AbstractArrayCollection

Array collection for carbonate deposition configuration.

# Fields
- `grid_carb::`[`CarbArray2DState`](@ref): `(7, toponum)` Carbonate grid array
    - ROW 1: x-coordinates of carbonate grid nodes.
    - ROW 2: Elevation at carbonate grid nodes.
    - ROW 3: Not used.
    - ROW 4: x-component of velocity at carbonate grid nodes.
    - ROW 5: y-component of velocity at carbonate grid nodes.
    - ROW 6: Not used.
    - ROW 7: Carbonate nucleation probability.

# Nested Dot Access
- `$(ROOT_NAME).grid_carb.array`

# Constructor
    Arrays(toponum::Int)::Arrays

Create a new Arrays collection with the given carbonate grid resolution.

# Arguments
- `toponum::Int`: Number of carbonate grid points
"""
mutable struct Arrays <: AbstractArrayCollection
    grid_carb::CarbArray2DState
end

function Arrays(toponum::Int)::Arrays
    units = [
        "1: m",
        "2: m",
        "3: m",
        "4: m/s",
        "5; m/s",
        "6: m",
        "7: m"
    ]
    descriptions = [
        "1: x-coordinates of carbonate grid nodes.",
        "2: Elevation at carbonate grid nodes.",
        "3: Not used.",
        "4: x-component of velocity at carbonate grid nodes.",
        "5: y-component of velocity at carbonate grid nodes.",
        "6: Not used.",
        "7: Carbonate nucleation probability."
    ]
    grid_carb = CarbArray2DState(
        zeros(Float64, 7, toponum),
        "grid_carb",
        units,
        descriptions
    )
    return Arrays(grid_carb)
end

end # module 