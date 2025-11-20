module TopoTypeManager

include("runge_kutta/RungeKutta.jl")
include("antidiffusion/Antidiffusion.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import ..Options: option_names

function topography_advect!(
    model::ModelData,
    inside_flags::Vector{Int8},
    ::Val{option_names.Antidiffusion}
)::Nothing
    @timeit_memit "Finished updating topography using velocity field (antidiffusion)" begin
        Antidiffusion.topography_advect!(model, inside_flags)
    end
    return nothing
end

function topography_advect!(
    model::ModelData,
    inside_flags::Vector{Int8},
    ::Val{option_names.RungeKuttaWithInterp}
)::Nothing
    @timeit_memit "Finished updating topography using velocity field (Runge-Kutta with interpolation)" begin
        RungeKutta.topography_advect!(model, inside_flags)
    end
    return nothing
end

end
