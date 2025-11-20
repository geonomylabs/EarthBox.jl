module PlotVelocityArraysManager

mutable struct PlotVelocityArrays
    gridx::Vector{Float64}
    gridy::Vector{Float64}
    vectorx::Matrix{Float64}
    vectory::Matrix{Float64}
    vectormag::Matrix{Float64}
    velocity_h5_datanames::Vector{String}

    function PlotVelocityArrays()
        new(
            zeros(Float64, 1),
            zeros(Float64, 1),
            zeros(Float64, 1, 1),
            zeros(Float64, 1, 1),
            zeros(Float64, 1, 1),
            String[]
        )
    end
end

function set_grid_arrays!(
    arrays::PlotVelocityArrays,
    gridx::Vector{Float64},
    gridy::Vector{Float64}
)::Nothing
    arrays.gridx = gridx
    arrays.gridy = gridy
    return nothing
end

function set_velocity_arrays!(
    arrays::PlotVelocityArrays,
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