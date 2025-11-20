module CompactionTools

""" Compact/decompact a layer using Athy's Law.

# Arguments
- `porosity_initial::Float64`: Initial porosity of the layer (fraction)
- `depth_decay_term::Float64`: Depth decay term (1/m)
- `top_initial::Vector{Float64}`: Initial top depth of the layer (m)
- `thickness_initial::Vector{Float64}`: Initial thickness of the layer (m)
- `top_new::Vector{Float64}`: New top depth of the layer (m)

# Returns
- `thickness_new::Vector{Float64}`: New thickness of the layer (m)
"""
function compact_or_decompact_layer(
    porosity_initial::Float64,
    depth_decay_term::Float64,
    top_initial::Vector{Float64},
    thickness_initial::Vector{Float64},
    top_new::Vector{Float64}
)::Vector{Float64}
    xnum = length(thickness_initial)
    thickness_new = Vector{Float64}(undef, xnum) #zeros(Float64, xnum)
    for i in 1:xnum
        thickness_new[i] = compact_or_decompact(
            porosity_initial,
            depth_decay_term,
            top_initial[i],
            top_initial[i] + thickness_initial[i],
            top_new[i]
        )
    end
    return thickness_new
end

function compact_or_decompact_layer(
    porosity_initial::Float64,
    depth_decay_term::Float64,
    top_initial::Float64,
    thickness_initial::Vector{Float64},
    top_new::Float64
)::Vector{Float64}
    xnum = length(thickness_initial)
    thickness_new = Vector{Float64}(undef, xnum) #zeros(Float64, xnum)
    for i in 1:xnum
        thickness_new[i] = compact_or_decompact(
            porosity_initial,
            depth_decay_term,
            top_initial,
            top_initial + thickness_initial[i],
            top_new
        )
    end
    return thickness_new
end


""" Compact or decompact a layer using Athy's Law.

# Arguments
- `porosity_initial::Float64`: Initial porosity of the layer (fraction)
- `depth_decay_term::Float64`: Depth decay term (1/m)
- `top_initial::Float64`: Initial top depth of the layer (m)
- `bottom_initial::Float64`: Initial bottom depth of the layer (m)
- `top_new::Float64`: New top depth of the layer (m)
- `nmax::Int`: Maximum number of iterations (default: 5)
- `tolerance::Float64`: Tolerance for convergence (default: 1e-6)

# Returns
- `thickness_new::Float64`: New thickness of the layer (m)
"""
function compact_or_decompact(
    porosity_initial::Float64,
    depth_decay_term::Float64,
    top_initial::Float64,
    bottom_initial::Float64,
    top_new::Float64;
    nmax::Int = 5,
    tolerance::Float64 = 1e-6
)::Float64
    thickness_initial = bottom_initial - top_initial
    thickness_old = thickness_initial
    bottom_new = top_new + thickness_old
    bottom_new_initial = bottom_new

    convergence_criterion = 1e32
    icount = 0

    while convergence_criterion > tolerance && icount < nmax
        exceeded_exponential_input_limit = check_exponential_input_limit(
            depth_decay_term, top_initial, bottom_initial, top_new, bottom_new)
        if !exceeded_exponential_input_limit
            exp1 = calculate_exponential_term_for_compaction(
                depth_decay_term, porosity_initial, top_initial, bottom_initial)
            exp2 = calculate_exponential_term_for_compaction(
                depth_decay_term, porosity_initial, top_new, bottom_new)
            thickness_new = thickness_initial - exp1 + exp2
            bottom_new = top_new + thickness_new
            convergence_criterion = abs(thickness_new - thickness_old)
            thickness_old = thickness_new
            icount += 1
        else
            convergence_criterion = tolerance
            bottom_new = bottom_new_initial
        end
    end

    thickness_new = bottom_new - top_new
    if thickness_new < 0.0
        thickness_new = 0.001
    end

    return thickness_new
end

""" Check if the exponential inputs are within the limits.

# Arguments
- `depth_decay_term::Float64`: Depth decay term (1/m)
- `top_initial::Float64`: Initial top depth (m)
- `bottom_initial::Float64`: Initial bottom depth (m)
- `top_new::Float64`: New top depth (m)
- `bottom_new::Float64`: New bottom depth (m)

# Returns
- `check::Bool`: Whether any input exceeds the limit
"""
function check_exponential_input_limit(
    depth_decay_term::Float64,
    top_initial::Float64,
    bottom_initial::Float64,
    top_new::Float64,
    bottom_new::Float64
)::Bool
    exponential_input_limit = 100.0
    check = false
    if abs(depth_decay_term * top_initial) > exponential_input_limit
        check = true
    end
    if abs(depth_decay_term * bottom_initial) > exponential_input_limit
        check = true
    end
    if abs(depth_decay_term * top_new) > exponential_input_limit
        check = true
    end
    if abs(depth_decay_term * bottom_new) > exponential_input_limit
        check = true
    end
    return check
end

""" Calculate the exponential term for compaction.

# Arguments
- `depth_decay_term::Float64`: Depth decay term (1/m)
- `porosity_initial::Float64`: Initial porosity (fraction)
- `top_depth::Float64`: Top depth (m)
- `bottom_depth::Float64`: Bottom depth (m)

# Returns
- `exponential_term::Float64`: Calculated exponential term
"""
function calculate_exponential_term_for_compaction(
    depth_decay_term::Float64,
    porosity_initial::Float64,
    top_depth::Float64,
    bottom_depth::Float64
)::Float64
    exponential_term = (
        porosity_initial / depth_decay_term * (
            exp(-depth_decay_term * top_depth) -
            exp(-depth_decay_term * bottom_depth)
        )
    )
    return exponential_term
end

""" Calculate final thickness of layer after burial.

# Arguments
- `porosity_initial::Float64`: Initial porosity (fraction)
- `depth_decay_term::Float64`: Depth decay term (1/m)
- `thickness_initial_x::Vector{Float64}`: Initial thicknesses (m)
- `total_thickness_final_x::Vector{Float64}`: Final total thicknesses (m)

# Returns
- `thickness_x::Vector{Float64}`: Final thicknesses (m)
"""
function calculate_final_thickness_after_burial(
    porosity_initial::Float64,
    depth_decay_term::Float64,
    thickness_initial_x::Vector{Float64},
    total_thickness_final_x::Vector{Float64}
)::Vector{Float64}
    xnum = length(thickness_initial_x)
    thickness_x = zeros(Float64, xnum)
    for i in 1:xnum
        thickness_x[i] = calculate_final_thickness(
            porosity_initial,
            depth_decay_term,
            thickness_initial_x[i],
            total_thickness_final_x[i]
        )
    end
    return thickness_x
end

""" Calculate final thickness of a layer after burial.

# Arguments
- `porosity_initial::Float64`: Initial porosity (fraction)
- `depth_decay_term::Float64`: Depth decay term (1/m)
- `thickness_initial::Float64`: Initial thickness (m)
- `total_thickness_final::Float64`: Final total thickness (m)

# Returns
- `thickness::Float64`: Final thickness (m)
"""
function calculate_final_thickness(
    porosity_initial::Float64,
    depth_decay_term::Float64,
    thickness_initial::Float64,
    total_thickness_final::Float64
)::Float64
    tolerance = 1e-5
    nmax = 10
    thickness = thickness_initial
    thickness_load_old = 0.0
    convergence_criterion = 1e32

    icount = 0
    while convergence_criterion > tolerance && icount <= nmax
        thickness_load_new = total_thickness_final - thickness
        thickness = compact_or_decompact(
            porosity_initial,
            depth_decay_term,
            thickness_load_old,
            thickness_load_old + thickness,
            thickness_load_new
        )
        convergence_criterion = abs(thickness_load_old - thickness_load_new)
        thickness_load_old = thickness_load_new
        icount += 1
    end
    return thickness
end

""" Print iteration information.

# Arguments
- `icount::Int`: Iteration count
- `thickness_load_old::Float64`: Old thickness load
- `thickness_load_new::Float64`: New thickness load
- `thickness::Float64`: Current thickness
- `convergence_criterion::Float64`: Convergence criterion
"""
function print_iteration_info(
    icount::Int,
    thickness_load_old::Float64,
    thickness_load_new::Float64,
    thickness::Float64,
    convergence_criterion::Float64
)::Nothing
    println(
        "Iteration ", icount,
        " thickness_load_old ", round(thickness_load_old, digits=5),
        " thickness_load_new ", round(thickness_load_new, digits=5),
        " thickness ", thickness,
        " convergence_criterion ", convergence_criterion
    )
end

end # module 