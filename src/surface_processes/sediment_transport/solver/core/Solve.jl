module Solve

import SparseArrays: SparseMatrixCSC, spzeros
import LinearAlgebra: lu
import EarthBox.MathTools: linear_interp_bisection
import EarthBox.DataStructures: SedimentTransportParameters
import ..BuildSys: build_sys_topo!
import ..BuildSys: build_sys_topo_tridiagonal!

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
    L::Matrix{Float64},
    R::Vector{Float64},
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

    build_sys_topo!(
        L, R,
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

""" Allocation-light tridiagonal variant of `solve_downhill_diffusion`.

Uses three preallocated diagonal vectors (`dl`, `d`, `du`) and a
preallocated solution buffer `S` instead of a dense (toponum × toponum)
`L::Matrix{Float64}` and `Ls = SparseMatrixCSC(L); lu(Ls) \\ R`. The
underlying mathematical system is identical (the legacy build is purely
tridiagonal), so the Thomas algorithm produces the same answer to within
floating-point rounding.

The Thomas algorithm is destructive on `d`, `du`, and `R`: forward
elimination overwrites them with the LU factors; back-substitution
writes the solution into `S`. The next call rebuilds `dl`/`d`/`du`/`R`
from scratch via `build_sys_topo_tridiagonal!`, so the destructive
overwrite is harmless.

# Returns
- `S`: same Vector{Float64} that was passed in, now containing the
    solution. Returning the same buffer (rather than allocating a new
    Vector{Float64}) keeps the per-call allocation at zero.
"""
function solve_downhill_diffusion_optimized(
    dl::Vector{Float64},
    d::Vector{Float64},
    du::Vector{Float64},
    R::Vector{Float64},
    S::Vector{Float64},
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

    build_sys_topo_tridiagonal!(
        dl, d, du, R,
        topo_gridx, topo_gridy,
        topo_grid_diffusivity,
        topo_grid_pelagic_sedimentation_rate,
        xmin_bgrid, xmax_bgrid,
        timestep,
        porosity_initial_pelagic, depth_decay_term_pelagic
    )

    thomas_solve!(dl, d, du, R, S)
    return S
end

""" Hand-rolled Thomas algorithm for tridiagonal systems.

Solves `T * x = R` where `T = Tridiagonal(dl, d, du)`, writing the
solution into preallocated buffer `S`. Forward sweep destructively
modifies `d` and `R` (they hold the running LU factors after the call);
back-substitution then fills `S`. `dl` is read but not modified.

Equivalent in floating-point to LAPACK `gtsv!` (which Julia's
`Tridiagonal \\ R` dispatches to) but written inline here so the per-call
allocation is zero — `Tridiagonal(dl,d,du) \\ R` allocates copies of all
inputs to preserve them.
"""
function thomas_solve!(
    dl::Vector{Float64},
    d::Vector{Float64},
    du::Vector{Float64},
    R::Vector{Float64},
    S::Vector{Float64}
)::Nothing
    n = length(d)
    @assert length(dl) == n - 1
    @assert length(du) == n - 1
    @assert length(R) == n
    @assert length(S) == n
    # Forward elimination: combine row i with row i-1 to zero out the
    # sub-diagonal at row i. d and R are mutated to hold the running LU
    # factor and modified RHS.
    @inbounds for i in 2:n
        m = dl[i-1] / d[i-1]
        d[i] -= m * du[i-1]
        R[i] -= m * R[i-1]
    end
    # Back substitution: solve the upper-bidiagonal system into S.
    @inbounds S[n] = R[n] / d[n]
    @inbounds for i in (n-1):-1:1
        S[i] = (R[i] - du[i] * S[i+1]) / d[i]
    end
    return nothing
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