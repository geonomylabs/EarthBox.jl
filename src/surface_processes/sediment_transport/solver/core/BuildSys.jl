"""
    BuildSys

Functions used to build system of equations for 1D sediment transport.
"""
module BuildSys

using EarthBox.Compaction.CompactionTools: compact_or_decompact

""" Build system of equations for downhill diffusion.

Equation for downhill diffusion:
    dYt/dt = d(K*dYt/dx)/dx + H

where:
    Yt = topography elevation (meters),
    K = diffusivity (m^2/s),
    H = background pelagic sedimentation rate (m/s),
    t = time (seconds),
    x = distance (meters).

1D grid used to discretize the equation:
    YtL-----YtC-----YtR
    KnL     KnC     KnR

FD representation for constant diffusivity (KnL = KnC = KnR = K):
    [YtC - YtC_old]/dt = K*[(YtR - YtC)/dx - (YtC - YtL)/dx]/dx + H

    YtC - K*dt*(YtL - 2*YtC + YtR)/dx^2 = YtC_old + H

    [1 + 2Kdt/dx^2]*YtC - K*dt/dx^2*YtL - K*dt/dx^2*YtR = YtC_old + H

FD representation for variable diffusivity:
    [YtC - YtC_old]/dt = [KR(YtR - YtC)/dx - KL(YtC - YtL)/dx]/dx

    YtC + (KL + KR)*dt/dx^2*YtC - KL*dt/dx^2*YtL - KR*dt/dx^2*YtR = YtC_old

    [1 + (KL + KR)*dt/dx^2]*YtC - KL*dt/dx^2*YtL - KR*dt/dx^2*YtR = YtC_old

    where KL = (KnL + KnC)/2 and KR = (KnC + KnR)/2.
"""
function build_sys_topo(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_grid_diffusivity::Vector{Float64},
    topo_grid_pelagic_sedimentation_rate::Vector{Float64},
    xmin_model_grid::Float64,
    xmax_model_grid::Float64,
    timestep::Float64,
    porosity_initial_pelagic::Float64,
    depth_decay_term_pelagic::Float64
)::Tuple{Matrix{Float64}, Vector{Float64}}
    toponum = length(topo_gridx)
    dx_topo = topo_gridx[2] - topo_gridx[1]
    L = zeros(toponum, toponum)
    R = zeros(toponum)
    
    apply_symmetry_boundary_condition_left_boundary(L, R, 1)
    
    for i in 2:toponum-1
        xtopo = topo_gridx[i]
        if inside_model_domain(xtopo, xmin_model_grid, xmax_model_grid)
            apply_internal_stencil_variable_diffusivity(
                L, R, topo_gridy, topo_grid_diffusivity,
                topo_grid_pelagic_sedimentation_rate, timestep, dx_topo, i,
                porosity_initial_pelagic, depth_decay_term_pelagic
            )
        else
            apply_symmetry_to_external_nodes(L, R, xtopo, xmin_model_grid, i)
        end
    end
    
    apply_symmetry_boundary_condition_right_boundary(L, R, toponum)
    return L, R
end

""" Check if topography node is inside model domain.

# Arguments
- `xtopo`: Topography node x-coordinate
- `xmin`: Minimum x-coordinate of model domain
- `xmax`: Maximum x-coordinate of model domain

# Returns
- `check`: True if topography node is inside model domain, False otherwise
"""
function inside_model_domain(
    xtopo::Float64,
    xmin::Float64,
    xmax::Float64
)::Bool
    return xmin <= xtopo <= xmax
end

""" Apply symmetry boundary conditions to external nodes.

# Arguments
- `L`: Matrix of coefficients
- `R`: Right-hand side of the equation
- `xtopo`: Topography node x-coordinate
- `xmin`: Minimum x-coordinate of model domain
- `i`: Index of the topography node
"""
function apply_symmetry_to_external_nodes(
    L::Matrix{Float64},
    R::Vector{Float64},
    xtopo::Float64,
    xmin::Float64,
    i::Int
)::Nothing
    if xtopo < xmin
        apply_symmetry_boundary_condition_left_boundary(L, R, i)
    else
        apply_symmetry_boundary_condition_right_boundary(L, R, i)
    end
    return nothing
end

""" Apply internal stencil.

The source term involves a decompacted pelagic sedimentation rate and is used to 
calculate a decompacted pelagic sediment thickness. However, the transport model 
works in terms of compacted sediment thickness. Therefore, the decompacted pelagic 
sediment thickness is converted to compacted sediment thickness using the input 
initial porosity and depth decay terms.

# Arguments
- `L`: Matrix of coefficients
- `R`: Right-hand side of the equation
- `topo_gridy`: Topography elevation (meters)
- `topo_grid_diffusivity`: Topography diffusion coefficient (m^2/s)
- `topo_grid_pelagic_sedimentation_rate`: Pelagic sedimentation rate (m/s)
- `timestep`: Time step (seconds)
- `dx_topo`: Topography grid spacing (meters)
- `i`: Index of the topography node
- `porosity_initial_pelagic`: Initial porosity of pelagic sediment at the 
    sediment-water interface (fraction)
- `depth_decay_term_pelagic`: Depth decay term of pelagic sediment (1/m)
"""
function apply_internal_stencil_variable_diffusivity(
    L::Matrix{Float64},
    R::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_grid_diffusivity::Vector{Float64},
    topo_grid_pelagic_sedimentation_rate::Vector{Float64},
    timestep::Float64,
    dx_topo::Float64,
    i::Int,
    porosity_initial_pelagic::Float64,
    depth_decay_term_pelagic::Float64
)::Nothing
    diffusivity_left_mid, diffusivity_right_mid = 
        calculate_average_midpoint_diffusivity(topo_grid_diffusivity, i)
    
    L[i, i-1] = -diffusivity_left_mid * timestep / dx_topo^2
    L[i, i] = 1.0 + (diffusivity_left_mid + diffusivity_right_mid) * 
        timestep / dx_topo^2
    L[i, i+1] = -diffusivity_right_mid * timestep / dx_topo^2
    
    compacted_thickness_pelagic = 
        calculate_the_compacted_thickness_for_pelagic_sediment(
            topo_grid_pelagic_sedimentation_rate[i],
            timestep, porosity_initial_pelagic, depth_decay_term_pelagic
        )
    R[i] = topo_gridy[i] - compacted_thickness_pelagic
    return nothing
end

function calculate_the_compacted_thickness_for_pelagic_sediment(
    sedimentation_rate::Float64,
    timestep::Float64,
    porosity_initial::Float64,
    depth_decay_term::Float64
)::Float64
    thickness = sedimentation_rate * timestep
    thickness_compacted = compact_or_decompact(
        porosity_initial, depth_decay_term, 0.0, thickness, 10_000.0)
    return thickness_compacted
end

""" Calculate average midpoint diffusivity.

# Arguments
- `topo_grid_diffusivity`: Topography diffusion coefficient (m^2/s)
- `i`: Index of the topography node

# Returns
- `diffusivity_left_mid`: Average diffusivity between left and center nodes
- `diffusivity_right_mid`: Average diffusivity between center and right nodes
"""
function calculate_average_midpoint_diffusivity(
    topo_grid_diffusivity::Vector{Float64},
    i::Int
)::Tuple{Float64, Float64}
    diffusivity_left_node = topo_grid_diffusivity[i-1]
    diffusivity_center_node = topo_grid_diffusivity[i]
    diffusivity_right_node = topo_grid_diffusivity[i+1]

    diffusivity_left_mid = 0.5 * (diffusivity_left_node + diffusivity_center_node)
    diffusivity_right_mid = 0.5 * (diffusivity_center_node + diffusivity_right_node)
    return diffusivity_left_mid, diffusivity_right_mid
end

""" Apply symmetry boundary conditions along left boundary.

This boundary condition involves zero flux across the left boundary.

# Arguments
- `L`: Matrix of coefficients
- `R`: Right-hand side of the equation
- `i`: Index of the topography node
"""
function apply_symmetry_boundary_condition_left_boundary(
    L::Matrix{Float64},
    R::Vector{Float64},
    i::Int
)::Nothing
    L[i, i] = 1
    L[i, i+1] = -1
    R[i] = 0
    return nothing
end

""" Apply symmetry boundary conditions along right boundary.

This boundary condition involves zero flux across the right boundary.

# Arguments
- `L`: Matrix of coefficients
- `R`: Right-hand side of the equation
- `i`: Index of the topography node
"""
function apply_symmetry_boundary_condition_right_boundary(
    L::Matrix{Float64},
    R::Vector{Float64},
    i::Int
)::Nothing
    L[i, i] = 1
    L[i, i-1] = -1
    R[i] = 0
    return nothing
end

end # module 