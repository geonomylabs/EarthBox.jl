module InterfaceDepthsContainer

""" Struct used to store interface depths.

Each depth array is defined at each basic grid x-coordinate and has xnum 
elements. All depths are in meters.
"""
mutable struct InterfaceDepths
    # Y-coordinates of the top of sticky layer (air and water)
    ytops_sticky::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the bottom of sticky layer (air and water)
    ybottoms_sticky::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the top of sticky layer
    ytops_sticky_air::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the bottom of sticky layer
    ybottoms_sticky_air::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the top of sticky water layer
    ytops_sticky_water::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the bottom of sticky water layer
    ybottoms_sticky_water::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the top of rock layer
    ytops_rock::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the bottom of rock layer
    ybottoms_rock::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the top of crust layer
    ytops_crust::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the bottom of crust layer
    ybottoms_crust::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the top of mantle lithosphere layer
    ytops_mantle_lith::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the bottom of mantle lithosphere layer
    ybottoms_mantle_lith::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the top of asthenosphere layer
    ytops_asthenosphere::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the bottom of asthenosphere layer
    ybottoms_asthenosphere::Union{Vector{Float64}, Nothing}
    # Y-Coordinate of topography
    ytopo::Union{Vector{Float64}, Nothing}
    # Y-Coordinate of Moho
    ymoho::Union{Vector{Float64}, Nothing}
    # Y-coordinates of smoothed topography
    ytopo_smooth::Union{Vector{Float64}, Nothing}
    # Y-Coordinate of base of lithosphere
    ylith_base::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the top of sediments
    ytops_sediments::Union{Vector{Float64}, Nothing}
    # Y-coordinates of the bottom of sediments
    ybottoms_sediments::Union{Vector{Float64}, Nothing}
end

# Constructor with default values
InterfaceDepths(; ytops_sticky=nothing, ybottoms_sticky=nothing,
    ytops_sticky_air=nothing, ybottoms_sticky_air=nothing,
    ytops_sticky_water=nothing, ybottoms_sticky_water=nothing,
    ytops_rock=nothing, ybottoms_rock=nothing,
    ytops_crust=nothing, ybottoms_crust=nothing,
    ytops_mantle_lith=nothing, ybottoms_mantle_lith=nothing,
    ytops_asthenosphere=nothing, ybottoms_asthenosphere=nothing,
    ytopo=nothing, ymoho=nothing, ytopo_smooth=nothing,
    ylith_base=nothing, ytops_sediments=nothing,
    ybottoms_sediments=nothing) = 
    InterfaceDepths(ytops_sticky, ybottoms_sticky, ytops_sticky_air,
        ybottoms_sticky_air, ytops_sticky_water, ybottoms_sticky_water,
        ytops_rock, ybottoms_rock, ytops_crust, ybottoms_crust,
        ytops_mantle_lith, ybottoms_mantle_lith, ytops_asthenosphere,
        ybottoms_asthenosphere, ytopo, ymoho, ytopo_smooth,
        ylith_base, ytops_sediments, ybottoms_sediments)

function calc_ytopo!(data::InterfaceDepths, xnum::Int)
    data.ytopo = calc_topography(xnum, data.ytops_rock, data.ybottoms_sticky)
end

function calc_ytopo_smooth!(data::InterfaceDepths, nsmooth::Int, xnum::Int)
    data.ytopo_smooth = smooth_interface(nsmooth, xnum, data.ytopo)
end

function calc_ymoho!(data::InterfaceDepths, xnum::Int)
    data.ymoho = calc_moho(xnum, data.ytops_mantle_lith, data.ybottoms_crust)
end

function calc_ylith_base!(data::InterfaceDepths, xnum::Int)
    data.ylith_base = calc_base_of_lithosphere(xnum, data.ytops_asthenosphere,
        data.ybottoms_mantle_lith)
end

function clean_crustal_depths!(data::InterfaceDepths, xnum::Int)
    clean_top_and_bottom_of_crust!(xnum, data.ytops_crust, data.ybottoms_crust,
        data.ytopo)
end

function clean_lithosphere_depths!(data::InterfaceDepths, xnum::Int)
    clean_top_and_bottom_of_mantle_lithosphere!(xnum, data.ytops_mantle_lith,
        data.ybottoms_mantle_lith, data.ybottoms_crust)
end

function clean_asthenosphere_depths!(data::InterfaceDepths, xnum::Int)
    clean_top_and_bottom_of_asthenosphere!(xnum, data.ytops_asthenosphere,
        data.ybottoms_asthenosphere, data.ymoho)
end

function calculate_water_depth(data::InterfaceDepths)
    ytopos_water_final = (data.ybottoms_sticky_air .+ data.ytops_sticky_water) ./ 2.0
    ybottoms_water_final = (data.ybottoms_sticky_water .+ data.ytops_rock) ./ 2.0
    return ybottoms_water_final .- ytopos_water_final
end

function print_depths(data::InterfaceDepths)
    println("ylith_base: ", data.ylith_base[3])
    println("ytops_mantle_lith, ybottoms_mantle_lith: ",
        data.ytops_mantle_lith[3], " ", data.ybottoms_mantle_lith[3])
    println("ytops_asthenosphere, ybottoms_mantle_lith: ",
        data.ytops_asthenosphere[3], " ", data.ybottoms_mantle_lith[3])
end

function calc_topography(
    xnum::Int,
    ytops_rock::Vector{Float64},
    ybottoms_sticky_water::Vector{Float64}
)::Vector{Float64}
    ytopo = zeros(Float64, xnum)
    for i in 1:xnum
        ytopo[i] = (ytops_rock[i] + ybottoms_sticky_water[i]) / 2.0
    end
    return ytopo
end

function smooth_topography(
    xnum::Int,
    ytopo::Vector{Float64}
)::Vector{Float64}
    ytopo_smooth = zeros(Float64, xnum)
    for i in 1:xnum
        if 1 < i < xnum
            ytopo_smooth[i] = (ytopo[i] + ytopo[i-1] + ytopo[i+1]) / 3.0
        elseif i == 1
            ytopo_smooth[i] = (ytopo[i] + ytopo[i+1]) / 2.0
        elseif i == xnum
            ytopo_smooth[i] = (ytopo[i] + ytopo[i-1]) / 2.0
        end
    end
    return ytopo_smooth
end

function smooth_interface(
    nsmooth::Int,
    xnum::Int,
    ytopo::Vector{Float64}
)::Vector{Float64}
    ytopo_smooth = zeros(Float64, xnum)
    for i in 1:xnum
        ytopo_sum = 0.0
        icount = 0
        for j in 0:(2*nsmooth)
            ii = i - nsmooth + j
            if 1 <= ii <= xnum
                ytopo_sum += ytopo[ii]
                icount += 1
            end
        end
        ytopo_smooth[i] = ytopo_sum / icount
    end
    return ytopo_smooth
end

function clean_top_and_bottom_of_crust!(
    xnum::Int,
    ytops_crust::Vector{Float64},
    ybottoms_crust::Vector{Float64},
    ytopo::Vector{Float64}
)::Nothing
    for i in 1:xnum
        if ytops_crust[i] == 0 && ybottoms_crust[i] == 0.0
            ytops_crust[i] = ytopo[i]
            ybottoms_crust[i] = ytopo[i]
        end
    end
end

function calc_moho(
    xnum::Int,
    ytops_mantle_lith::Vector{Float64},
    ybottoms_crust::Vector{Float64}
)::Vector{Float64}
    ymoho = zeros(Float64, xnum)
    for i in 1:xnum
        ymoho[i] = (ytops_mantle_lith[i] + ybottoms_crust[i]) / 2.0
    end
    return ymoho
end

function calc_base_of_lithosphere(
    xnum::Int,
    ytops_asthenosphere::Vector{Float64},
    ybottoms_mantle_lith::Vector{Float64}
)::Vector{Float64}
    ylith_base = zeros(Float64, xnum)
    for i in 1:xnum
        ylith_base[i] = (ytops_asthenosphere[i] + ybottoms_mantle_lith[i]) / 2.0
    end
    return ylith_base
end

function clean_top_and_bottom_of_mantle_lithosphere!(
    xnum::Int,
    ytops_mantle_lith::Vector{Float64},
    ybottoms_mantle_lith::Vector{Float64},
    ybottoms_crust::Vector{Float64}
)::Nothing
    for i in 1:xnum
        if ytops_mantle_lith[i] == 0 && ybottoms_mantle_lith[i] == 0.0
            ytops_mantle_lith[i] = ybottoms_crust[i]
            ybottoms_mantle_lith[i] = ybottoms_crust[i]
        end
    end
end

function clean_top_and_bottom_of_asthenosphere!(
    xnum::Int,
    ytops_asthenosphere::Vector{Float64},
    ybottoms_asthenosphere::Vector{Float64},
    ymoho::Vector{Float64}
)::Nothing
    for i in 1:xnum
        if ytops_asthenosphere[i] == 0 && ybottoms_asthenosphere[i] == 0.0
            ytops_asthenosphere[i] = ymoho[i]
            ybottoms_asthenosphere[i] = ymoho[i]
        end
    end
end

end # module 