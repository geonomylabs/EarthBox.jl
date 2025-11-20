module PlotScalarArraysManager

import ..DataNames: ScalarH5Datanames

mutable struct PlotScalarArrays
    gridx::Vector{Float64}
    gridy::Vector{Float64}
    scalar::Matrix{Float64}
    scalar_data_names::ScalarH5Datanames
end

function PlotScalarArrays()
    return PlotScalarArrays(
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1, 1),
        ScalarH5Datanames()
    )
end

function set_grid_arrays!(
    arrays::PlotScalarArrays, 
    gridx::Vector{Float64}, 
    gridy::Vector{Float64}
)::Nothing
    arrays.gridx = gridx
    arrays.gridy = gridy
    return nothing
end

function set_scalar_arrays!(
    arrays::PlotScalarArrays, 
    scalar::Matrix{Float64}
)::Nothing
    arrays.scalar = copy(scalar)
    return nothing
end

end # module 