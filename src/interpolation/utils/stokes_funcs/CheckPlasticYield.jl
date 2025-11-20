module CheckPlasticYield

"""
    plastic_yielding(
        i::Int, 
        j::Int, 
        plasticn::Matrix{Float64}, 
        plastics::Matrix{Float64}, 
        plastic_yield::Matrix{Float64}, 
        plas_tol::Float64
    )::Bool

Check for plastic yielding on basic grid nodes.

# Arguments
- `i::Int`: First index
- `j::Int`: Second index
- `plasticn::Matrix{Float64}`: Array of plastic values
- `plastics::Matrix{Float64}`: Array of plastic values
- `plastic_yield::Matrix{Float64}`: Array of plastic yield values
- `plas_tol::Float64`: Plastic tolerance threshold

# Returns
- `Bool`: True if plastic yielding is detected, false otherwise
"""
@inline function plastic_yielding(
    i::Int, 
    j::Int, 
    plasticn::Matrix{Float64}, 
    plastics::Matrix{Float64}, 
    plastic_yield::Matrix{Float64}, 
    plas_tol::Float64
)::Bool
    check = false
    
    if plasticn[i,j] > plas_tol
        check = true
    end
    
    if (plastics[i,j] > plas_tol ||
        plastics[i,j+1] > plas_tol ||
        plastics[i+1,j] > plas_tol ||
        plastics[i+1,j+1] > plas_tol)
        check = true
    end
    
    if (plastic_yield[i,j] > plas_tol ||
        plastic_yield[i,j+1] > plas_tol ||
        plastic_yield[i+1,j] > plas_tol ||
        plastic_yield[i+1,j+1] > plas_tol)
        check = true
    end
    
    return check
end

end # module 