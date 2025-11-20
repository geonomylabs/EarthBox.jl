module PlotVectorArraysManager

mutable struct PlotVectorArrays
    gridx::Vector{Float64}
    gridy::Vector{Float64}
    vectorx::Matrix{Float64}
    vectory::Matrix{Float64}
    vectormag::Matrix{Float64}
    velocity_h5_datanames::Vector{String}
end

function PlotVectorArrays()
    return PlotVectorArrays(
        zeros(Float64, 1),
        zeros(Float64, 1),
        zeros(Float64, 1, 1),
        zeros(Float64, 1, 1),
        zeros(Float64, 1, 1),
        ["velx_cmyr", "vely_cmyr", "velmag_cmyr"]
    )
end

function set_grid_arrays!(
    arrays::PlotVectorArrays, 
    gridx::Vector{Float64}, 
    gridy::Vector{Float64}
)::Nothing
    arrays.gridx = gridx
    arrays.gridy = gridy
    return nothing
end

function set_vector_arrays!(
    arrays::PlotVectorArrays,
    vectorx::Matrix{Float64},
    vectory::Matrix{Float64},
    vectormag::Matrix{Float64}
)::Nothing
    arrays.vectorx = copy(vectorx)
    arrays.vectory = copy(vectory)
    arrays.vectormag = copy(vectormag)
    return nothing
end

end # module 