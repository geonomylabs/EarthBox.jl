module WeightFuncs3d

@inline function upper_left_back_node_weight(
    dxULB::Float64,
    dyULB::Float64,
    dzULB::Float64,
)::Float64
    return (1.0-dxULB) * (1.0-dyULB) * (1.0-dzULB)
end

@inline function upper_left_front_node_weight(
    dxULB::Float64,
    dyULB::Float64,
    dzULB::Float64,
)::Float64
    return (1.0-dxULB) * (1.0-dyULB) * dzULB
end

@inline function lower_left_back_node_weight(
    dxULB::Float64,
    dyULB::Float64,
    dzULB::Float64,
)::Float64
    return (1.0-dxULB) * dyULB * (1.0-dzULB)
end

@inline function lower_left_front_node_weight(
    dxULB::Float64,
    dyULB::Float64,
    dzULB::Float64,
)::Float64
    return (1.0-dxULB) * dyULB * dzULB
end

@inline function upper_right_back_node_weight(
    dxULB::Float64,
    dyULB::Float64,
    dzULB::Float64,
)::Float64
    return dxULB * (1.0-dyULB) * (1.0-dzULB)
end

@inline function upper_right_front_node_weight(
    dxULB::Float64,
    dyULB::Float64,
    dzULB::Float64,
)::Float64
    return dxULB * (1.0-dyULB) * dzULB
end

@inline function lower_right_back_node_weight(
    dxULB::Float64,
    dyULB::Float64,
    dzULB::Float64,
)::Float64
    return dxULB * dyULB * (1.0-dzULB)
end

@inline function lower_right_front_node_weight(
    dxULB::Float64,
    dyULB::Float64,
    dzULB::Float64,
)::Float64
    return dxULB * dyULB * dzULB
end

end
