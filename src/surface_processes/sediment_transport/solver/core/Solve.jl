module Solve

import SparseArrays: SparseMatrixCSC, spzeros
import LinearAlgebra: lu
import EarthBox.MathTools: linear_interp_bisection
import EarthBox.DataStructures: SedimentTransportParameters
import ..BuildSys: build_sys_topo

""" Solve downhill diffusion model.

# Arguments
- `topo_gridx`: X-coordinates of topography grid (meters)
- `topo_gridy`: Y-coordinates of topography grid (meters)
- `topo_grid_diffusivity`: Diffusivity grid (m²/yr)
- `topo_grid_pelagic_sedimentation_rate`: Pelagic sedimentation rate grid (m/yr)
- `basic_grid_x_dimensions`: Tuple of (xmin, xmax) for basic grid (meters)
- `timestep`: Time step (years)
- `sediment_transport_parameters`: Sediment transport parameters

# Returns
- `S`: Solution vector
"""
function solve_downhill_diffusion(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_grid_diffusivity::Vector{Float64},
    topo_grid_pelagic_sedimentation_rate::Vector{Float64},
    basic_grid_x_dimensions::Tuple{Float64, Float64},
    timestep::Float64,
    sediment_transport_parameters::SedimentTransportParameters
)::Vector{Float64}
    porosity_initial_pelagic = sediment_transport_parameters.porosity_initial
    depth_decay_term_pelagic = sediment_transport_parameters.depth_decay_term
    
    xmin_bgrid = basic_grid_x_dimensions[1]
    xmax_bgrid = basic_grid_x_dimensions[2]
    
    L, R = build_sys_topo(
        topo_gridx, topo_gridy,
        topo_grid_diffusivity,
        topo_grid_pelagic_sedimentation_rate,
        xmin_bgrid, xmax_bgrid,
        timestep,
        porosity_initial_pelagic, depth_decay_term_pelagic
    )
    
    Ls = SparseMatrixCSC(L)
    S = lu(Ls) \ R
    return S
end

""" Update elevation for vertical advection and diffusion.

# Arguments
- `topo_gridx_new`: New topography grid x-coordinates (meters)
- `topo_gridy_new`: New topography grid y-coordinates (meters)
- `gridt`: Multi-dimensional topography array. Y-coordinate (meters) of
    topography nodes are stored in gridt[2, :].

# Updated Arrays
- `gridt`: Multi-dimensional topography array is updated with interpolated y-coordinates
"""
function interpolate_new_topography_to_topography_array(
    topo_gridx_new::Vector{Float64},
    topo_gridy_new::Vector{Float64},
    gridt::Matrix{Float64}
)::Nothing
    toponum = length(topo_gridy_new)
    for i in 1:toponum
        x_location = gridt[1, i]
        y_interp = linear_interp_bisection(
            topo_gridx_new, topo_gridy_new, x_location
        )
        gridt[2, i] = y_interp
    end
    
    return nothing
end

end # module 