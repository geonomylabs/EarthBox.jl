module PlotGridArraysManager

mutable struct PlotGridArrays
    gridx::Vector{Float64}
    gridy::Vector{Float64}
    gs::Matrix{Float64}
    vx::Matrix{Float64}
    vy::Matrix{Float64}
    vmag::Matrix{Float64}
    velocity_h5_datanames::Vector{String}
end

function PlotGridArrays()
    return PlotGridArrays(
        zeros(Float64, 1), # gridx
        zeros(Float64, 1), # gridy
        zeros(Float64, 1, 1), # gs
        zeros(Float64, 1, 1), # vx
        zeros(Float64, 1, 1), # vy
        zeros(Float64, 1, 1), # vmag
        ["velx_cmyr", "vely_cmyr", "velmag_cmyr"]
    )
end

function set_grid_arrays!(
    arrays::PlotGridArrays, 
    gridx::Vector{Float64}, 
    gridy::Vector{Float64}
)::Nothing
    arrays.gridx = gridx
    arrays.gridy = gridy
    return nothing
end

function set_scalar_arrays!(
    arrays::PlotGridArrays, 
    dataname::String, 
    ebarray2d::Matrix{Float64}
)::Nothing
    if dataname == "velx_cmyr"
        set_vx_array!(arrays, copy(ebarray2d))
    elseif dataname == "vely_cmyr"
        set_vy_array!(arrays, copy(ebarray2d))
    elseif dataname == "velmag_cmyr"
        set_vmag_array!(arrays, copy(ebarray2d))
    else
        set_grid_scalar_array!(arrays, copy(ebarray2d))
    end
    return nothing
end

function set_grid_scalar_array!(
    arrays::PlotGridArrays, 
    scalar::Matrix{Float64}
)::Nothing
    arrays.gs = scalar
    return nothing
end

function set_vx_array!(
    arrays::PlotGridArrays, 
    vx::Matrix{Float64}
)::Nothing
    arrays.vx = vx
    return nothing
end

function set_vy_array!(
    arrays::PlotGridArrays, 
    vy::Matrix{Float64}
)::Nothing
    arrays.vy = vy
    return nothing
end

function set_vmag_array!(
    arrays::PlotGridArrays, 
    vmag::Matrix{Float64}
)::Nothing
    arrays.vmag = vmag
    return nothing
end

end # module 