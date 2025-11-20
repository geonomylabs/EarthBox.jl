module SolverManager

include("utils/GlobalIndices.jl")
include("norms/StokesNorms.jl")
include("rhs/StokesContinuityRhs.jl")
include("rhs/StokesViscoelasticTerms.jl")
include("residuals/StokesResiduals.jl")
include("build_system/StokesBuildManager.jl")
include("system_solver/SystemSolver.jl")

end # module

