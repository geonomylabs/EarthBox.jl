module Profiles

import ..Reader
using Plots

""" Get a profile of scalar values at x-index icol.
"""
function get_yprofile_scalar(
    input_path::String,
    base_name::String,
    ioutput::Int,
    icol::Int
)::Tuple{Vector{Float64}, Vector{Float64}, Float64}
    gridx, gridy, vals, tMyr = Reader.read_scalar_grid(
        input_path, base_name, ioutput)
    ynum = length(gridy)
    profile = zeros(ynum)
    for j in 1:ynum
        profile[j] = vals[j, icol]
    end
    return gridy, profile, tMyr
end

""" Get a profile of scalar values at y-index irow.
"""
function get_xprofile_scalar(
    input_path::String,
    base_name::String,
    ioutput::Int, 
    irow::Int
)::Tuple{Vector{Float64}, Vector{Float64}, Float64}
    gridx, gridy, vals, tMyr = Reader.read_scalar_grid(
        input_path, base_name, ioutput)
    xnum = length(gridx)
    profile = zeros(xnum)
    for j in 1:xnum
        profile[j] = vals[irow, j]
    end
    return gridx, profile, tMyr
end

""" Get a profile of vx values at y-index irow.
"""
function get_xprofile_vx(
    input_path::String,
    base_name::String,
    ioutput::Int,
    irow::Int
)::Tuple{Vector{Float64}, Vector{Float64}, Float64}
    gridx, gridy, vx, vy, tmyr = Reader.read_vector_grid(
        input_path, base_name, ioutput)
    xnum = length(gridx)
    profile = zeros(xnum)
    for j in 1:xnum
        profile[j] = vx[irow, j]
    end
    return gridx, profile, tmyr
end

""" Get a profile of vy values at y-index irow.
"""
function get_xprofile_vy(
    input_path::String,
    base_name::String,
    ioutput::Int,
    irow::Int
)::Tuple{Vector{Float64}, Vector{Float64}, Float64}
    gridx, gridy, vx, vy, tmyr = Reader.read_vector_grid(
        input_path, base_name, ioutput)
    xnum = length(gridx)
    profile = zeros(xnum)
    for j in 1:xnum
        profile[j] = vy[irow, j]
    end
    return gridx, profile, tmyr
end

""" Plot temperature profile at specified x-index.
"""
function plot_temperature_profile(input_path::String, ioutput::Int, jcheck::Int)
    println(">> Reading temperature file ")
    base_name = "TempC"

    gridx, gridy, tk1_C, tMyr = Reader.read_scalar_grid(
        input_path, base_name, ioutput)
    
    xnum = length(gridx)
    ynum = length(gridy)
    snum = length(tk1_C)

    T_min = minimum(tk1_C)
    T_max = maximum(tk1_C)

    xmin = gridx[1]
    xmax = gridx[end]
    ymin = gridy[1]
    ymax = gridy[end]

    println("xmin, xmax, ymin, ymax, T_min, T_max : ",
            xmin, " ", xmax, " ", ymin, " ", ymax, " ", T_min, " ", T_max)
    println("xnum, ynum, snum : ", xnum, " ", ynum, " ", snum)

    println(" >> Extract temperature profile")

    temp_profile = zeros(ynum)
    j = jcheck
    xm = gridx[j] * 1000.0

    println(" ****** Extracting information at j, xm : ", j, " ", xm)

    for i in 1:ynum
        temp_profile[i] = tk1_C[i, j]
    end

    println(" >> Plotting temperature profile ")

    p = plot(gridy, temp_profile, seriestype=:scatter, color=:red)
    
    base_name = "TempC_profile"
    plot_name = base_name * "_" * string(tMyr) * "_" * string(ioutput) * ".png"
    savefig(p, plot_name)
end

""" Plot conductivity profile at specified x-index.
"""
function plot_conductivity_profile(input_path::String, ioutput::Int, jcheck::Int)
    println(">> Reading conductivity file ")
    base_name = "therm_cond"
    gridx, gridy, k, tMyr = Reader.read_scalar_grid(
        input_path, base_name, ioutput)

    xnum = length(gridx)
    ynum = length(gridy)
    snum = length(k)

    k_min = minimum(k)
    k_max = maximum(k)

    xmin = gridx[1]
    xmax = gridx[end]
    ymin = gridy[1]
    ymax = gridy[end]

    println("xmin, xmax, ymin, ymax, k_min, k_max : ",
            xmin, " ", xmax, " ", ymin, " ", ymax, " ", k_min, " ", k_max)
    println("xnum, ynum, snum : ", xnum, " ", ynum, " ", snum)

    println(" >> Extract conductivity profile")

    k_profile = zeros(ynum)
    j = jcheck
    xm = gridx[j] * 1000.0

    println(" ****** Extracting information at j, xm : ", j, " ", xm)

    for i in 1:ynum
        k_profile[i] = k[i, j]
    end

    println(" >> Plotting conductivity profile ")

    p = plot(gridy, k_profile, seriestype=:scatter, color=:red)
    
    base_name = "k_profile"
    plot_name = base_name * "_" * string(tMyr) * "_" * string(ioutput) * ".png"
    savefig(p, plot_name)
end

end  # end of module