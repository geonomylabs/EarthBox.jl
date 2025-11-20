module YieldCheck

@inline function check_for_yielding(
    no_yielding_in_mobile_wall::Bool,
    mat_id::Int16,
    matid_sticky_air::Int16,
    matid_mobile_wall::Int16,
    matid_plate_extension::Int16,
    stress_yield::Float64,
    stress_invariant_forecast::Float64
)::Tuple{Bool, Int}
    if no_yielding_in_mobile_wall
        yielding, yielding_flag = check_for_yielding_sandbox_mobile_wall(
            mat_id, matid_sticky_air, matid_mobile_wall, matid_plate_extension,
            stress_yield, stress_invariant_forecast
        )
    else
        yielding, yielding_flag = check_for_yielding_for_general_case(
            stress_yield, stress_invariant_forecast
        )
    end
    return yielding, yielding_flag
end

@inline function check_for_yielding_for_general_case(
    stress_yield::Float64,
    stress_invariant::Float64
)::Tuple{Bool, Int}
    yielding = false
    yielding_flag = 0
    if stress_invariant > stress_yield
        yielding = true
        yielding_flag = 1
    end
    return yielding, yielding_flag
end

""" Check for yielding for sandbox model with mobile wall. 

There is also an optional plate extension associated with the mobile wall for 
sandbox extension experiments.
"""
@inline function check_for_yielding_sandbox_mobile_wall(
    mat_id::Int16,
    matid_sticky_air::Int16,
    matid_mobile_wall::Int16,
    matid_plate_extension::Int16,
    stress_yield::Float64,
    stress_invariant::Float64
)::Tuple{Bool, Int}
    yielding = false
    yielding_flag = 0
    # Correcting rock properties for yielding
    # Do not allow yielding in mobile wall or sticky air
    if stress_invariant > stress_yield &&
       mat_id != matid_sticky_air &&
       mat_id != matid_mobile_wall &&
       mat_id != matid_plate_extension
        yielding = true
        yielding_flag = 1
    end
    return yielding, yielding_flag
end

end # module YieldCheck 