""" Test the calculation of the gravity anomaly of a cell.

The analytical solution for an infinite beam from Telford et al. (1990) is
compared with the solution from Turcotte and Schubert (2014) for a finite
prism.

"""
module TestCellGravity

import EarthBox.Gravity.GravityCell: calc_cell_gravity_for_xarray, calculate_gravity_anomaly_prism
import EarthBox.ConversionFuncs: meters_per_second_squared_to_milligal
import Plots

""" Box geometry for the gravity cell.
"""
mutable struct BoxGeometry
    x1_box::Float64
    x2_box::Float64
    y1_box::Float64
    y2_box::Float64
    z1_box::Float64
    z2_box::Float64

    function BoxGeometry(
        x1_box::Float64, x2_box::Float64,
        y1_box::Float64, y2_box::Float64,
        z1_box::Float64, z2_box::Float64
    )
        new(x1_box, x2_box, y1_box, y2_box, z1_box, z2_box)
    end
end

""" Calculate the width of the box.
"""
function width(box::BoxGeometry)::Float64
    return box.x2_box - box.x1_box
end

""" Calculate the depth of the box.
"""
function depth(box::BoxGeometry)::Float64
    return box.y2_box - box.y1_box
end

""" Calculate the height of the box.
"""
function height(box::BoxGeometry)::Float64
    return box.z2_box - box.z1_box
end

""" Calculate the x-coordinate of the center of the box.
"""
function xmid(box::BoxGeometry)::Float64
    return (box.x1_box + box.x2_box) / 2.0
end

mutable struct ObserverGeometry
    xmin_observer::Float64
    xmax_observer::Float64
    xspacing_observer::Float64
    xcoors_observer::Vector{Float64}
    ycoors_observer::Vector{Float64}
    zcoors_observer::Vector{Float64}

    function ObserverGeometry(
        xmin_observer::Float64, 
        xmax_observer::Float64,
        xspacing_observer::Float64,
        ycoor_observer::Float64, 
        zcoor_observer::Float64
    )
        xcoors = calc_xcoors(xmin_observer, xmax_observer, xspacing_observer)
        new(xmin_observer, xmax_observer, xspacing_observer, xcoors,
            [ycoor_observer], [zcoor_observer])
    end
end

""" Calculate the number of x-observers.
"""
function xnum(obs::ObserverGeometry)::Int
    xnum_val = floor(Int, (obs.xmax_observer - obs.xmin_observer) / 
        obs.xspacing_observer) + 1
    return xnum_val
end

""" Calculate the x-coordinates of the observers.
"""
function calc_xcoors(
    xmin_observer::Float64, 
    xmax_observer::Float64, 
    xspacing_observer::Float64
)::Vector{Float64}
    xnum_val = floor(Int, (xmax_observer - xmin_observer) / 
        xspacing_observer) + 1
    println("xnum = $xnum_val")
    coors = zeros(Float64, xnum_val)
    for i in 1:xnum_val
        xcoor_observer = xmin_observer + Float64(i - 1) * xspacing_observer
        coors[i] = xcoor_observer
    end
    return coors
end

""" Gravity anomaly of 3d prism.
"""
mutable struct GravityPrism
    gravity_constant::Float64
    box_geometry::BoxGeometry
    observer_geometry::ObserverGeometry
    delta_density::Float64
    gravity_anomaly_mgal::Vector{Float64}
    gravity_anomaly_nondim::Vector{Float64}

    function GravityPrism(
        box_geometry::BoxGeometry,
        observer_geometry::ObserverGeometry,
        delta_density::Float64
    )
        new(6.6732e-11, box_geometry, observer_geometry, delta_density,
            Float64[], Float64[])
    end
end

""" Calculate the gravity anomaly.
"""
function calculate_gravity(prism::GravityPrism)::Nothing
    gravity_constant = prism.gravity_constant
    delta_density = prism.delta_density

    x1_box = prism.box_geometry.x1_box
    x2_box = prism.box_geometry.x2_box
    y1_box = prism.box_geometry.y1_box
    y2_box = prism.box_geometry.y2_box
    z1_box = prism.box_geometry.z1_box
    z2_box = prism.box_geometry.z2_box

    xnum_val = xnum(prism.observer_geometry)
    for i in 1:xnum_val
        xcoor_observer = prism.observer_geometry.xcoors_observer[i]
        dx1 = x1_box - xcoor_observer
        dx2 = x2_box - xcoor_observer

        ycoor_observer = prism.observer_geometry.ycoors_observer[1]
        dy1 = y1_box - ycoor_observer
        dy2 = y2_box - ycoor_observer

        zcoor_observer = prism.observer_geometry.zcoors_observer[1]
        dz1 = z1_box - zcoor_observer
        dz2 = z2_box - zcoor_observer

        delta_gravity = calculate_gravity_anomaly_prism(
            gravity_constant,
            delta_density,
            dx1, dx2,
            dy1, dy2,
            dz1, dz2
        )

        delta_gravity_mgal = meters_per_second_squared_to_milligal(
            delta_gravity)

        push!(prism.gravity_anomaly_mgal, delta_gravity_mgal)
        push!(prism.gravity_anomaly_nondim, 
            delta_gravity / 2.0 / gravity_constant / delta_density)
    end
end

""" Plot the gravity anomaly in mGal.
"""
function plot_gravity_mgal(prism::GravityPrism)::Nothing
    p = Plots.plot(
        prism.observer_geometry.xcoors_observer, 
        prism.gravity_anomaly_mgal)
    Plots.xlabel!(p, "x_observer (m)")
    Plots.ylabel!(p, "gravity anomaly(mgl)")
    Plots.title!(p, "Test gravprism function")
    Plots.display(p)
end

""" Plot the gravity anomaly in non-dimensional units.
"""
function plot_gravity_nondim(prism::GravityPrism)::Nothing
    p = Plots.plot(
        prism.observer_geometry.xcoors_observer, 
        prism.gravity_anomaly_nondim)
    Plots.xlabel!(p, "x_observer (m)")
    Plots.ylabel!(p, "gravity anomaly / 2 / G / delta_density")
    Plots.title!(p, "Test gravprism function")
    Plots.display(p)
end

""" Calculate gravity anomaly of an infinite prism.
"""
mutable struct GravityInfinitePrism
    gravity_constant::Float64
    box_geometry::BoxGeometry
    observer_geometry::ObserverGeometry
    delta_density::Float64
    gravity_anomaly_nondim::Union{Vector{Float64}, Nothing}

    function GravityInfinitePrism(
        box_geometry::BoxGeometry,
        observer_geometry::ObserverGeometry,
        delta_density::Float64
    )
        new(6.6732e-11, box_geometry, observer_geometry, delta_density, 
            nothing)
    end
end

""" Calculate the gravity anomaly.
"""
function calculate_gravity(prism::GravityInfinitePrism)::Nothing
    cell_width = width(prism.box_geometry)
    cell_height = height(prism.box_geometry)
    z_datum = 0.0
    x_upper_left = prism.box_geometry.x1_box
    z_upper_left = prism.box_geometry.z1_box
    
    grav_anomaly_x = calc_cell_gravity_for_xarray(
        prism.gravity_constant,
        prism.delta_density,
        x_upper_left, z_upper_left, z_datum,
        cell_height, cell_width,
        prism.observer_geometry.xcoors_observer,
    )
    prism.gravity_anomaly_nondim = (
        grav_anomaly_x / (2.0 * prism.gravity_constant * prism.delta_density))
    return nothing
end

""" Plot the gravity anomaly in non-dimensional units.
"""
function plot_nondim(prism::GravityInfinitePrism)::Nothing
    p = Plots.plot(
        prism.observer_geometry.xcoors_observer,
        prism.gravity_anomaly_nondim
    )
    Plots.xlabel!(p, "x_observe (m)")
    Plots.ylabel!(p, "gravity_anomaly / 2 / G / delta_density")
    Plots.title!(p, "Gravity anomaly of a cell below a datum")
    Plots.display(p)
end

""" Test parameters for the gravity cell.
"""
mutable struct TestParameters
    box_geometry::Union{BoxGeometry, Nothing}
    observer_geometry::Union{ObserverGeometry, Nothing}
    delta_density::Union{Float64, Nothing}

    function TestParameters()
        new(nothing, nothing, nothing)
    end
end

""" Set the parameters from Turcotte and Schubert (2014).

Chapter 12, Problem 6.

"""
function set_parameters_from_turcotte14_ch12p6!(params::TestParameters)::Nothing
    params.delta_density = 400.0 # kg/m^3

    # observer location
    _xmin_observer = -20_000.0 # meters
    _xmax_observer = 20_000.0 # meters
    _xspacing_observer = 500.0 # meters
    _ycoor_observer = 0.0
    _zcoor_observer = 0.0

    params.observer_geometry = ObserverGeometry(
        _xmin_observer,
        _xmax_observer,
        _xspacing_observer,
        _ycoor_observer,
        _zcoor_observer
    )

    # prism location
    xwidth_box = 10_000.0 # meters
    _x1_box = -xwidth_box / 2.0
    _x2_box = xwidth_box / 2.0
    _y1_box = -50_000.0 # meters
    _y2_box = 50_000.0 # meters
    _z1_box = 1_000.0 # meters
    _z2_box = 1_500.0 # meters

    params.box_geometry = BoxGeometry(
        x1_box=_x1_box, x2_box=_x2_box,
        y1_box=_y1_box, y2_box=_y2_box,
        z1_box=_z1_box, z2_box=_z2_box
    )
end

""" Set parameters from Telford et al. (1990).
"""
function set_parameters_from_telfordf2p31!(params::TestParameters)::Nothing
    params.delta_density = 1.0 # kg/m^3

    # observer location
    _xmin_observer = -3.0 # meters
    _xmax_observer = 3.0 # meters
    _xspacing_observer = 0.1 # meters
    _ycoor_observer = 0.0
    _zcoor_observer = 0.0

    params.observer_geometry = ObserverGeometry(
        _xmin_observer,
        _xmax_observer,
        _xspacing_observer,
        _ycoor_observer,
        _zcoor_observer
    )

    # prism location
    xwidth_box = 1.0 # meters
    _x1_box = -xwidth_box / 2.0
    _x2_box = xwidth_box / 2.0
    _y1_box = -100.0 # meters
    _y2_box = 100.0 # meters
    _z1_box = 1.0 / 3.0 # meters
    _z2_box = 4.0 / 3.0 # meters

    params.box_geometry = BoxGeometry(
        _x1_box, _x2_box,
        _y1_box, _y2_box,
        _z1_box, _z2_box
    )
    return nothing
end

""" Compare infinite prism solution to finite prism solution.
"""
function run_test()::Nothing
    parameters = TestParameters()
    set_parameters_from_telfordf2p31!(parameters)

    box_geometry = parameters.box_geometry
    observer_geometry = parameters.observer_geometry
    delta_density = parameters.delta_density

    gravity_infinite_prism = GravityInfinitePrism(
        box_geometry,
        observer_geometry,
        delta_density
    )
    calculate_gravity(gravity_infinite_prism)
    gravity_anomaly_telford = gravity_infinite_prism.gravity_anomaly_nondim

    gravity_prism = GravityPrism(
        box_geometry,
        observer_geometry,
        delta_density
    )
    calculate_gravity(gravity_prism)
    gravity_anomaly_turcotte = gravity_prism.gravity_anomaly_nondim

    xcoors_observer = observer_geometry.xcoors_observer
    p = Plots.plot(
        xcoors_observer, gravity_anomaly_turcotte,
        label="Turcotte", color=:blue, linewidth=4
    )
    Plots.plot!(
        p, xcoors_observer, gravity_anomaly_telford,
        label="Telford", color=:red, linestyle=:dash
    )
    Plots.xlabel!(p, "x_observer (m)")
    Plots.ylabel!(p, "gravity anomaly / 2 / G / delta_density")
    Plots.title!(p, "Gravity anomaly of a cell below a datum")
    #Plots.legend!(p)
    Plots.display(p)
end

end # module
