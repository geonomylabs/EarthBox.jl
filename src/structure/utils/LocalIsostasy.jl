module LocalIsostasy

""" Calculate y coordinate of local isostatic topography.

Returns
-------
ytopo_local::Vector{Float64}
    Y-coordinates of local isostatic topography relative to a left-hand
    reference location.
"""
function calc_local_isostatic_topography_for_2dgrid(
    xindex_reference::Int,
    x_reference::Float64,
    y_bottom::Float64,
    xnum::Int,
    ytopo::Vector{Float64},
    ymoho::Vector{Float64},
    ylith_base::Vector{Float64},
    rho_sticky_x::Vector{Float64},
    rho_crust_x::Vector{Float64},
    rho_lithosphere_x::Vector{Float64},
    rho_asthenosphere_x::Vector{Float64}
)::Vector{Float64}
    gravitational_acceleration = 9.8 # m/s/s

    # Index of reference used to calculate local isostatic force balance
    # This index defines the reference lithosphere column
    rho_mantle_ref = rho_asthenosphere_x[xindex_reference]
    rho_sticky_ref = rho_sticky_x[xindex_reference]
    pressure_ref = calculate_pressure_of_column_for_index(
        xindex_reference,
        y_bottom, ytopo, ymoho, ylith_base,
        rho_sticky_x, rho_crust_x,
        rho_lithosphere_x, rho_asthenosphere_x,
        gravitational_acceleration
    )

    print_info(
        xindex_reference, x_reference,
        "Reference Lithosphere (Local Isostatic Topography Using Markers):",
        ytopo, ymoho, ylith_base, y_bottom,
        rho_crust_x, rho_lithosphere_x, rho_asthenosphere_x, rho_sticky_x,
        pressure_ref
    )

    ytopo_local = calc_local_isostatic_topo(
        xnum, y_bottom, ytopo, ymoho, ylith_base,
        rho_sticky_x, rho_crust_x,
        rho_lithosphere_x, rho_asthenosphere_x,
        pressure_ref, rho_sticky_ref, rho_mantle_ref,
        gravitational_acceleration
    )

    return ytopo_local
end

""" Calculate local topography.
"""
function calc_local_isostatic_topo(
    xnum::Int,
    y_bottom::Float64,
    ytopo::Vector{Float64},
    ymoho::Vector{Float64},
    ylith_base::Vector{Float64},
    rho_sticky_x::Vector{Float64},
    rho_crust_x::Vector{Float64},
    rho_lithosphere_x::Vector{Float64},
    rho_asthenosphere_x::Vector{Float64},
    pressure_ref::Float64,
    rho_sticky_ref::Float64,
    rho_mantle_ref::Float64,
    gravitational_acceleration::Float64
)::Vector{Float64}
    ytopo_local = zeros(Float64, xnum)
    for i in 1:xnum
        if i > 1
            pressure_column = calculate_pressure_of_column_for_index(
                i, y_bottom, ytopo, ymoho, ylith_base,
                rho_sticky_x, rho_crust_x,
                rho_lithosphere_x, rho_asthenosphere_x,
                gravitational_acceleration
            )
            y_correction = calculate_topo_correction(
                pressure_ref,
                pressure_column,
                rho_sticky_ref,
                rho_mantle_ref,
                gravitational_acceleration
            )
        else
            y_correction = 0.0
        end
        ytopo_local[i] = ytopo[i] + y_correction
    end
    return ytopo_local
end

""" Calculate pressure of column for index.
"""
function calculate_pressure_of_column_for_index(
    ix::Int,
    y_bottom::Float64,
    ytopo::Vector{Float64},
    ymoho::Vector{Float64},
    ylith_base::Vector{Float64},
    rho_sticky_x::Vector{Float64},
    rho_crust_x::Vector{Float64},
    rho_lithosphere_x::Vector{Float64},
    rho_asthenosphere_x::Vector{Float64},
    gravitational_acceleration::Float64
)::Float64
    # Get reference thicknesses of layers in column
    thickness_sticky = ytopo[ix]
    thickness_crust = ymoho[ix] - ytopo[ix]
    thickness_mantle_lithosphere = ylith_base[ix] - ymoho[ix]
    thickness_asthenosphere = y_bottom - ylith_base[ix]
    
    # Get average densities of layers in column
    rho_sticky = rho_sticky_x[ix]
    rho_crust = rho_crust_x[ix]
    rho_mantle_lithosphere = rho_lithosphere_x[ix]
    rho_asthenosphere = rho_asthenosphere_x[ix]

    pressure = calculate_pressure_of_column(
        max(thickness_sticky, 0),
        max(thickness_crust, 0),
        max(thickness_mantle_lithosphere, 0),
        max(thickness_asthenosphere, 0),
        rho_sticky,
        rho_crust,
        rho_mantle_lithosphere,
        rho_asthenosphere,
        gravitational_acceleration
    )
    return pressure
end

""" Calculate pressure of column.

Inputs
------
thickness_sticky::Float64
    Thickness of sticky air layer in m.
thickness_crust::Float64
    Thickness of crust in m.
thickness_mantle_lithosphere::Float64
    Thickness of mantle lithosphere in m.
thickness_asthenosphere::Float64
    Thickness of asthenosphere in m.
rho_sticky::Float64
    Density of sticky air layer in kg/m^3.
rho_crust::Float64
    Density of crust in kg/m^3.
rho_mantle_lithosphere::Float64
    Density of mantle lithosphere in kg/m^3.
rho_asthenosphere::Float64
    Density of asthenosphere in kg/m^3.

Returns
-------
pressure::Float64
    Pressure of column in Pa.
"""
function calculate_pressure_of_column(
    thickness_sticky::Float64,
    thickness_crust::Float64,
    thickness_mantle_lithosphere::Float64,
    thickness_asthenosphere::Float64,
    rho_sticky::Float64,
    rho_crust::Float64,
    rho_mantle_lithosphere::Float64,
    rho_asthenosphere::Float64,
    gravitational_acceleration::Float64
)::Float64
    pressure = (
        rho_sticky * thickness_sticky +
        rho_crust * thickness_crust +
        rho_mantle_lithosphere * thickness_mantle_lithosphere +
        rho_asthenosphere * thickness_asthenosphere
    ) * gravitational_acceleration
    return pressure
end

""" Calculate topography correction.

y_correction is the isostatic topography correction necessary to
balance the column with the reference
"""
function calculate_topo_correction(
    pressure_ref::Float64,
    pressure_column::Float64,
    rho_water::Float64,
    rho_mantle::Float64,
    gravitational_acceleration::Float64
)::Float64
    y_correction = (
        (pressure_ref - pressure_column) /
        (rho_water - rho_mantle) /
        gravitational_acceleration
    )
    return y_correction
end

function print_info(
    i::Int,
    x_reference::Float64,
    msg::String,
    ytopo::Vector{Float64},
    ymoho::Vector{Float64},
    ylith_base::Vector{Float64},
    y_bottom::Float64,
    rho_crust_x::Vector{Float64},
    rho_lithosphere_x::Vector{Float64},
    rho_asthenosphere_x::Vector{Float64},
    rho_sticky_x::Vector{Float64},
    basal_pressure::Float64
)::Nothing
    println(msg)
    println("x-index: $i")
    println("x-coordinate: $(x_reference/1000.0)")
    println("ytopo (km): $(ytopo[i]/1000.0)")
    println("ymoho (km): $(ymoho[i]/1000.0)")
    println("ylith (km): $(ylith_base[i]/1000.0)")
    println("y_bottom (km): $(y_bottom/1000.0)")
    println("sticky air + water thickness (km): $(ytopo[i]/1000.0)")
    println("crustal thickness (km): $((ymoho[i] - ytopo[i])/1000.0)")
    println("mantle lith. thick (km): $((ylith_base[i] - ymoho[i])/1000.0)")
    println("asthenos. thick (km): $((y_bottom - ylith_base[i])/1000.0)")
    println("rho_crust (kg/m^3): $(rho_crust_x[i])")
    println("rho_lithosphere (kg/m^3): $(rho_lithosphere_x[i])")
    println("rho_asthenosphere (kg/m^3): $(rho_asthenosphere_x[i])")
    println("rho_sticky (kg/m^3): $(rho_sticky_x[i])")
    println("basal pressure (GPa): $(basal_pressure/1e9)")
    return nothing
end

end # module 