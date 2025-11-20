module CheckYield

@inline function check_yield_is_applicable(
    ntimestep::Int,
    cohesion::Float64,
    sine_friction_angle::Float64
)::Bool
    if ntimestep > 0 && (cohesion > 0.0 || sine_friction_angle > 0.0)
        return true
    end
    return false
end

end # module CheckYield 