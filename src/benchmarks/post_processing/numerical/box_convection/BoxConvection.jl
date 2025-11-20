module BoxConvection

""" Perform post-processing calculations for convection in a box benchmark.

# Returns
- `vrms::Float64`: Root mean square of velocity in m/s
- `temp_avg::Float64`: Average temperature in Kelvins
- `nusselt_number::Float64`: Nusselt Number
"""
function box_convection_quantities(input_dict::Dict{String,Any})
    temperature_top = input_dict["temperature_top"]
    temperature_bottom = input_dict["temperature_bottom"]
    xnum = input_dict["xnum"]
    ynum = input_dict["ynum"]
    xsize = input_dict["xsize"]
    ysize = input_dict["ysize"]
    vx1 = input_dict["vx1"]
    vy1 = input_dict["vy1"]
    tk1 = input_dict["tk1"]
    xstp_b = input_dict["xstp_b"]
    ystp_b = input_dict["ystp_b"]
    density = input_dict["density"]
    heat_capacity = input_dict["heat_capacity"]
    conductivity = input_dict["conductivity"]

    vrms = calculate_nondimensional_rms_velocity(
        xnum, ynum, xsize, ysize,
        conductivity, heat_capacity, density,
        vx1, vy1, xstp_b, ystp_b
    )

    temp_avg = calculate_average_temperature(
        xnum, ynum, xsize, ysize,
        tk1, xstp_b, ystp_b
    )

    nusselt_number = calculate_nusselt_number(
        xnum, xsize, ysize,
        temperature_top, temperature_bottom,
        tk1, xstp_b, ystp_b
    )

    return vrms, temp_avg, nusselt_number
end

function calculate_average_temperature(
    xnum::Int,
    ynum::Int,
    xsize::Float64,
    ysize::Float64,
    tk1::Matrix{Float64},
    xstp_b::Vector{Float64},
    ystp_b::Vector{Float64}
)::Float64
    temp_avg = 0.0
    for j in 1:(xnum-1)
        for i in 1:(ynum-1)
            temp_avg += (
                tk1[i,j] + tk1[i+1,j] + tk1[i,j+1] + tk1[i+1,j+1]
            ) / 4.0 * xstp_b[j] * ystp_b[i]
        end
    end
    return temp_avg / xsize / ysize
end

function calculate_nondimensional_rms_velocity(
    xnum::Int,
    ynum::Int,
    xsize::Float64,
    ysize::Float64,
    conductivity::Float64,
    heat_capacity::Float64,
    density::Float64,
    vx1::Matrix{Float64},
    vy1::Matrix{Float64},
    xstp_b::Vector{Float64},
    ystp_b::Vector{Float64}
)::Float64
    vrms = 0.0
    for j in 1:(xnum-1)
        for i in 1:(ynum-1)
            vrms += (
                vx1[i+1,j]^2 + vx1[i+1,j+1]^2 + vy1[i,j+1]^2 + vy1[i+1,j+1]^2
            ) / 2.0 * xstp_b[j] * ystp_b[i]
        end
    end
    return sqrt(vrms / xsize / ysize) * ysize / (conductivity / heat_capacity / density)
end

function calculate_nusselt_number(
    xnum::Int,
    xsize::Float64,
    ysize::Float64,
    temperature_top::Float64,
    temperature_bottom::Float64,
    tk1::Matrix{Float64},
    xstp_b::Vector{Float64},
    ystp_b::Vector{Float64}
)::Float64
    nnus1 = 0.0
    for j in 1:(xnum-1)
        dtdy1 = (
            tk1[1,j] - tk1[2,j] + tk1[1,j+1] - tk1[2,j+1]
        ) / 2.0 / ystp_b[1]
        nnus1 -= dtdy1 * xstp_b[j]
    end
    return ysize * nnus1 / ((temperature_bottom - temperature_top) * xsize)
end

end # module 