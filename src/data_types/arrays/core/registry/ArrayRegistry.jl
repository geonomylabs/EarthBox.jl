module ArrayRegistry

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.Arrays.ArrayTypes.RhsHeatArray1D: RhsHeatArray1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.Arrays.ArrayTypes.Array1DInt: Array1DIntState
import EarthBox.Arrays.ArrayTypes.SolutionArray1D: SolutionArray1DState
import EarthBox.Arrays.ArrayTypes.RhsStokesArray1D: RhsStokesArray1DState
import EarthBox.Arrays.ArrayTypes.BcArrayFloat: BcArrayFloatState
import EarthBox.Arrays.ArrayTypes.InternalBcArrayInt: InternalBcArrayIntState
import EarthBox.Arrays.ArrayTypes.InternalBcArrayFloat: InternalBcArrayFloatState
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray2D
import EarthBox.Arrays.ArrayMacros: @arrays

struct ArrayData
    name::String
    units::String
    type::Union{Type{<:AbstractEarthBoxArray1D}, Type{<:AbstractEarthBoxArray2D}}
    grid_type::String
    description::String
end

# --- Section files ---------------------------------------------------------
# Each file under sections/ defines exactly one get_<subsystem>_arrays()
# function. The order below mirrors the merge order in get_eb_arrays() below.
include("sections/grid.jl")
include("sections/heat_equation.jl")
include("sections/interpolation.jl")
include("sections/markers.jl")
include("sections/melting.jl")
include("sections/stokes_continuity.jl")
include("sections/bcs.jl")

function get_eb_arrays()::NamedTuple
    return merge_named_tuples(
        get_grid_arrays(),
        get_heat_equation_arrays(),
        get_interpolation_arrays(),
        get_markers_arrays(),
        get_melting_arrays(),
        get_stokes_continuity_arrays(),
        get_bcs_arrays()
    )
end

function merge_named_tuples(tuples...)
    result = tuples[1]
    ntuples = length(tuples)
    for i in 2:ntuples
        result = merge(result, tuples[i])
    end
    return result
end

end # module
