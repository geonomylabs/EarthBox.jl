module Ticks

function get_colorbar_ticks_for_composition_plot(
    n_bin::Int
)::Vector{Float64}
    ticks = Float64[]
    for i in 1:n_bin
        push!(ticks, Float64(i))
    end
    return ticks
end

function get_colorbar_ticks(
    smax::Float64,
    nlevels::Int
)::Vector{Float64}
    dval = smax / nlevels
    ticks = Float64[]
    for i in 0:nlevels
        val = dval * Float64(i)
        push!(ticks, val)
    end
    return ticks
end

function get_colorbar_ticks_general(
    smin::Float64,
    smax::Float64,
    contour_interval::Float64
)::Vector{Float64}
    tick_interval = contour_interval * 2.0
    nlevels = calculate_number_of_levels(smin, smax, tick_interval)
    if nlevels > 50
        tick_interval = tick_interval * 4.0
        nlevels = calculate_number_of_levels(smin, smax, tick_interval)
    end
    ticks = Float64[]
    for i in 0:nlevels
        val = smin + tick_interval * Float64(i)
        push!(ticks, val)
    end
    return ticks
end

function calculate_number_of_levels(
    scalar_minimum::Float64,
    scalar_maximum::Float64,
    tick_interval::Float64
)::Int
    nlevels = floor(Int, (scalar_maximum - scalar_minimum) / tick_interval)
    return nlevels
end

function get_colorbar_ticks_age(
    scalar_minimum::Float64,
    scalar_maximum::Float64,
    contour_interval::Float64
)::Vector{Float64}
    nlevels = floor(Int, (scalar_maximum - scalar_minimum) / contour_interval)
    ticks = Float64[]
    for i in 0:nlevels
        val = scalar_minimum + contour_interval * Float64(i)
        push!(ticks, val)
    end
    return ticks
end

end # module