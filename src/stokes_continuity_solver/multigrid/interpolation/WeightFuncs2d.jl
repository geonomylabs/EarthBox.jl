module WeightFuncs2d

@inline function upper_left_node_weight(
    dxUL::Float64,
    dyUL::Float64,
)::Float64
    return (1.0-dxUL) * (1.0-dyUL)
end

@inline function lower_left_node_weight(
    dxUL::Float64,
    dyUL::Float64,
)::Float64
    return (1.0-dxUL) * dyUL
end

@inline function upper_right_node_weight(
    dxUL::Float64,
    dyUL::Float64,
)::Float64
    return dxUL * (1.0-dyUL)
end

@inline function lower_right_node_weight(
    dxUL::Float64,
    dyUL::Float64,
)::Float64
    return dxUL * dyUL
end

end
