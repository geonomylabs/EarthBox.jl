module ArrayTypes

include("types/Array1DInt.jl")
include("types/Array1DFloat.jl")
include("types/GridArray1D.jl")
include("types/MarkerArrayFloat1D.jl")
include("types/MarkerArrayInt1D.jl")
include("types/RhsStokesArray1D.jl")
include("types/RhsHeatArray1D.jl")
include("types/ScalarArray2D.jl")
include("types/ScalarArray3D.jl")
include("types/SolutionArray1D.jl")
include("types/BcArrayFloat.jl")
include("types/InternalBcArrayFloat.jl")
include("types/InternalBcArrayInt.jl")
include("types/TopoArray2D.jl")
include("types/CarbArray2D.jl")
include("types/MaterialArrayInt2D.jl")
include("types/MaterialArrayFloat2D.jl")
include("types/MaterialArrayFloat1D.jl")
include("types/MaterialArrayInt1D.jl")

end # module ArrayTypes