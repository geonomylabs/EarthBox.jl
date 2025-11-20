module SmoothingGroup

import EarthBox.Arrays.ArrayTypes.Array1DInt: Array1DIntState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

mutable struct Smoothing <: AbstractArrayGroup
    smoothing_iterations::Array1DIntState
end

function Smoothing(smoothing_iterations_on_finest_level::Int64, max_levels::Int64)::Smoothing
    smoothing_iterations = zeros(Int64, max_levels)
    set_smoothing_iterations(smoothing_iterations, smoothing_iterations_on_finest_level)
    return Smoothing(
        Array1DIntState(
            smoothing_iterations, "smoothing_iterations", "None",
            "Number of smoothing iterations used for Gauss-Seidel solver for each resolution level"
        )
    )
end

function set_smoothing_iterations(
    smoothing_iterations::Vector{Int64},
    smoothing_iterations_on_finest_level::Int64,
)::Nothing
    smoothing_iterations[1] = smoothing_iterations_on_finest_level
    max_levels = length(smoothing_iterations)
    for i = 2:max_levels
        smoothing_iterations[i] = smoothing_iterations[i-1] * 2
    end
    return nothing
end

end