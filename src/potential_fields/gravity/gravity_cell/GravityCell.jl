module GravityCell

function calc_cell_gravity_for_xarray(
        grav_constant::Float64,
        delta_density::Float64,
        x_upper_left::Float64,
        y_upper_left::Float64,
        y_datum::Float64,
        cell_height::Float64,
        cell_width::Float64,
        xarray::Vector{Float64}
)::Vector{Float64}
    nx = length(xarray)
    grav_anomaly_x = zeros(Float64, nx)
    for i in 1:nx
        x_observer = xarray[i]
        grav_anomaly_x[i] = calculate_gravity_anomaly_of_cell(
            grav_constant,
            delta_density,
            x_upper_left, y_upper_left,
            x_observer, y_datum,
            cell_height, cell_width
        )
    end
    return grav_anomaly_x
end

""" Calculate the gravity anomaly of a cell.

The cell has a density contrast `delta_density` relative to the reference medium. 
The upper left x-coordinate of the cell is `x_upper_left` and the sub-datum depth 
of the cell is `depth_upper_left`. The observer is located at x-coordinate 
`x_observer`.

# Arguments
- `grav_constant`: Gravitational constant (m³/kg/s²)
- `delta_density`: Density contrast of the cell relative to the reference medium (kg/m³)
- `x_upper_left`: x-coordinate (m) of the upper left corner of the cell
- `y_upper_left`: y-coordinate (m) of the upper left corner of the cell where y is 
  positive downward into the earth
- `x_observer`: x-coordinate of the observer (m)
- `y_datum`: y-coordinate of the datum (m)
- `cell_height`: Height of the cell (m)
- `cell_width`: Width of the cell (m)
- `beta`: Angle between the vertical and the side of the cell (degrees)

# Returns
- `Float64`: Gravity anomaly of the cell (mGal)
"""
function calculate_gravity_anomaly_of_cell(
        grav_constant::Float64,
        delta_density::Float64,
        x_upper_left::Float64,
        y_upper_left::Float64,
        x_observer::Float64,
        y_datum::Float64,
        cell_height::Float64,
        cell_width::Float64;
        beta::Float64 = 0.0
)::Float64
    (
        radius1, radius2, radius3, radius4, theta1, theta2, theta3, theta4
    ) = calculate_corner_radii_to_observer_and_angles_relative_to_vertical(
        x_upper_left, y_upper_left,
        x_observer, y_datum,
        cell_height, cell_width
    )
    x_relative_to_upper_left = x_observer - x_upper_left
    y_upper_left_below_datum = y_upper_left - y_datum
    y_lower_left_below_datum = y_upper_left_below_datum + cell_height

    beta_rad = beta * π / 180.0
    factor = calculate_gravity_factor(
        beta_rad, x_relative_to_upper_left,
        y_upper_left_below_datum, y_lower_left_below_datum, cell_width,
        radius1, radius2, radius3, radius4,
        theta1, theta2, theta3, theta4
    )

    gravity_anomaly = 2.0 * grav_constant * delta_density * factor

    return gravity_anomaly
end

function calculate_gravity_factor(
        beta::Float64,
        x_relative_to_upper_left::Float64,
        y_upper_left_below_datum::Float64,
        y_lower_left_below_datum::Float64,
        cell_width::Float64,
        radius1::Float64,
        radius2::Float64,
        radius3::Float64,
        radius4::Float64,
        theta1::Float64,
        theta2::Float64,
        theta3::Float64,
        theta4::Float64
)::Float64
    factor = (
        y_lower_left_below_datum * (theta2 - theta4)
        - y_upper_left_below_datum * (theta1 - theta3)
        + sin(beta) * cos(beta) * (
            x_relative_to_upper_left * (theta2 + theta3 - theta4 - theta1)
            + cell_width * (theta4 - theta3)
        )
        + cos(beta) * cos(beta) * (
            x_relative_to_upper_left * log(radius2 * radius3 / (radius4 * radius1))
            + cell_width * log(radius4 / radius3)
        )
    )
    return factor
end

function calculate_gravity_factor_when_beta_equal_zero(
        x_relative_to_upper_left::Float64,
        y_upper_left_below_datum::Float64,
        y_lower_left_below_datum::Float64,
        cell_width::Float64,
        radius1::Float64,
        radius2::Float64,
        radius3::Float64,
        radius4::Float64,
        theta1::Float64,
        theta2::Float64,
        theta3::Float64,
        theta4::Float64
)::Float64
    factor = (
        y_lower_left_below_datum * (theta2 - theta4)
        - y_upper_left_below_datum * (theta1 - theta3)
        + x_relative_to_upper_left * log(radius2 * radius3 / (radius4 * radius1))
        + cell_width * log(radius4 / radius3)
    )
    return factor
end

""" Calculate the corner radii and angles.

The radii are the distances from the observer to the corners of the cell.
The angles are the angles between the vertical and the line connecting the
observer to the corners of the cell.
"""
function calculate_corner_radii_to_observer_and_angles_relative_to_vertical(
        x_upper_left::Float64,
        y_upper_left::Float64,
        x_observer::Float64,
        y_datum::Float64,
        cell_height::Float64,
        cell_width::Float64
    )::Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}
    x_relative_to_upper_left = x_observer - x_upper_left
    y_upper_left_below_datum = y_upper_left - y_datum

    radius1, theta1 = corner_upper_left(
        x_relative_to_upper_left, y_upper_left_below_datum)

    radius2, theta2 = corner_lower_left(
        x_relative_to_upper_left, y_upper_left_below_datum, cell_height)

    radius3, theta3 = corner_upper_right(
        x_relative_to_upper_left, y_upper_left_below_datum, cell_width)

    radius4, theta4 = corner_lower_right(
        x_relative_to_upper_left, y_upper_left_below_datum,
        cell_height, cell_width)

    return radius1, radius2, radius3, radius4, theta1, theta2, theta3, theta4
end

""" Calculate the radius and angle of the upper left corner of the cell.
"""
function corner_upper_left(
        x_relative_to_upper_left::Float64,
        y_upper_left_below_datum::Float64
    )::Tuple{Float64, Float64}
    theta1 = calculate_angle_between_radius_and_vertical(
        x_relative_to_upper_left,
        y_upper_left_below_datum
    )
    radius1 = calc_corner_observer_radius(
        x_relative_to_upper_left,
        y_upper_left_below_datum
    )
    return radius1, theta1
end

""" Calculate the radius and angle of the lower left corner of the cell.
"""
function corner_lower_left(
        x_relative_to_upper_left::Float64,
        y_upper_left_below_datum::Float64,
        cell_height::Float64
    )::Tuple{Float64, Float64}
    theta2 = calculate_angle_between_radius_and_vertical(
        x_relative_to_upper_left,
        y_upper_left_below_datum + cell_height
    )
    radius2 = calc_corner_observer_radius(
        x_relative_to_upper_left,
        y_upper_left_below_datum + cell_height
    )
    return radius2, theta2
end

""" Calculate the radius and angle of the upper right corner of the cell.
"""
function corner_upper_right(
        x_relative_to_upper_left::Float64,
        y_upper_left_below_datum::Float64,
        cell_width::Float64
    )::Tuple{Float64, Float64}
    theta3 = calculate_angle_between_radius_and_vertical(
        (x_relative_to_upper_left - cell_width),
        y_upper_left_below_datum
    )
    radius3 = calc_corner_observer_radius(
        x_relative_to_upper_left - cell_width,
        y_upper_left_below_datum
    )
    return radius3, theta3
end

""" Calculate the radius and angle of the lower right corner of the cell.
"""
function corner_lower_right(
        x_relative_to_upper_left::Float64,
        y_upper_left_below_datum::Float64,
        cell_height::Float64,
        cell_width::Float64
    )::Tuple{Float64, Float64}
    theta4 = calculate_angle_between_radius_and_vertical(
        (x_relative_to_upper_left - cell_width),
        y_upper_left_below_datum + cell_height
    )
    radius4 = calc_corner_observer_radius(
        x_relative_to_upper_left - cell_width,
        y_upper_left_below_datum + cell_height
    )
    return radius4, theta4
end

function calculate_angle_between_radius_and_vertical(opposite::Float64, adjacent::Float64)::Float64
    if adjacent == 0.0
        theta = π / 2.0
    else
        theta = atan(opposite / adjacent)
    end
    return theta
end

function calc_corner_observer_radius(x_corner::Float64, y_corner::Float64)::Float64
    return sqrt(x_corner^2 + y_corner^2)
end

function print_info(
    x_observer::Float64, x_relative_to_upper_left::Float64,
    cell_width::Float64, grav_constant::Float64,
    delta_density::Float64,
    radius1::Float64, radius2::Float64, radius3::Float64, radius4::Float64,
    theta1::Float64, theta2::Float64, theta3::Float64, theta4::Float64,
    gravity_anomaly::Float64
)::Nothing
    println("x_observer: ", x_observer)
    println("x_relative_to_upper_left: ", x_relative_to_upper_left)
    println("x_relative_to_upper_left/cell_width: ", x_relative_to_upper_left/cell_width)
    println("radius1: ", radius1)
    println("radius2: ", radius2)
    println("radius3: ", radius3)
    println("radius4: ", radius4)
    println("theta1 (deg) ", radians_to_degrees(theta1))
    println("theta2 (deg) ", radians_to_degrees(theta2))
    println("theta3 (deg) ", radians_to_degrees(theta3))
    println("theta4 (deg) ", radians_to_degrees(theta4))
    println("gravity_anomaly ", gravity_anomaly)
    println("g/2/G/delta_density ", gravity_anomaly/2/grav_constant/delta_density)
end

function radians_to_degrees(theta::Float64)::Float64
    return theta * 180.0 / π
end

""" Calculate the gravity anomaly of a prism.

The geometry of the prism is defined relative to an observation point Pobs
(See Figure 1).

This code is a modified version from Turcotte and Schubert (2014), chapter 12.6.

               ----> +x
             /|
            / |
          +y  +z
                    o Pobs (xo, yo, zo)----------------------
                    |                                ^      ^
                    |     <----->                    |      |
                    |       dx1                      | dz1  |
                    |     <------------------->      |      |
                    |               dx2              V      |
                    |   ^       p1------------p2     -------|
                    |  /       / |            /|            | dz2
                    | /dy1    /  |           / |            |
                    |/______ /   |          /  |            |
                    /       /    |         /   |            V
                   /dy2    /    p5--------/---p6     -------
                  /       /    /         /    /
                 V       /    /         /    /
                        p3------------p4    /
                        |   /          |   /
                        |  /           |  /
                        | /            | /
                        |/             |/
                        p7------------p8
            Figure 1: Prism geometry and observation point.

Observer-node vectors
---------------------
Pobs --> p1: (dx1, dy1, dz1)
Pobs --> p2: (dx2, dy1, dz1)

Pobs --> p3: (dx1, dy2, dz1)
Pobs --> p4: (dx2, dy2, dz1)

Pobs --> p5: (dx1, dy1, dz2)
Pobs --> p6: (dx2, dy1, dz2)

Pobs --> p7: (dx1, dy2, dz2)
Pobs --> p8: (dx2, dy2, dz2)

# Arguments
- `gravity_constant`: Gravitational constant (Nm²/kg²)
- `delta_density`: Density contrast of the prism relative to the reference medium (kg/m³)
- `dx1`: x-coordinate of the left edge of the prism relative to the observation point (m)
- `dx2`: x-coordinate of the right edge of the prism relative to the observation point (m)
- `dy1`: y-coordinate of the forward of the prism relative to the observation point (m)
- `dy2`: y-coordinate of the back of the prism relative to the observation point (m)
- `dz1`: z-coordinate of the top of the prism relative to the observation point (m)
- `dz2`: z-coordinate of the bottom of the prism relative to the observation point (m)

# Returns
- `Float64`: Gravity anomaly of the prism (m/s²)
"""
function calculate_gravity_anomaly_prism(
        gravity_constant::Float64,
        delta_density::Float64,
        dx1::Float64,
        dx2::Float64,
        dy1::Float64,
        dy2::Float64,
        dz1::Float64,
        dz2::Float64
)::Float64
    # Work on the top points of the prism
    vector_mag_p1 = calculate_vector_magnitude(dx1, dy1, dz1)
    vector_mag_p2 = calculate_vector_magnitude(dx2, dy1, dz1)
    vector_mag_p3 = calculate_vector_magnitude(dx1, dy2, dz1)
    vector_mag_p4 = calculate_vector_magnitude(dx2, dy2, dz1)
    # Work on the bottom points of the prism
    vector_mag_p5 = calculate_vector_magnitude(dx1, dy1, dz2)
    vector_mag_p6 = calculate_vector_magnitude(dx2, dy1, dz2)
    vector_mag_p7 = calculate_vector_magnitude(dx1, dy2, dz2)
    vector_mag_p8 = calculate_vector_magnitude(dx2, dy2, dz2)
    # Calculate the gravity terms for top points of the prism
    gravity_term_p1 = calculate_gravity_term(
        dx1, dy1, dz1, vector_mag_p1, sign_factor=-1.0)
    gravity_term_p2 = calculate_gravity_term(
        dx2, dy1, dz1, vector_mag_p2, sign_factor=1.0)
    gravity_term_p3 = calculate_gravity_term(
        dx1, dy2, dz1, vector_mag_p3, sign_factor=1.0)
    gravity_term_p4 = calculate_gravity_term(
        dx2, dy2, dz1, vector_mag_p4, sign_factor=-1.0)
    # Calculate the gravity terms for bottom points of the prism
    gravity_term_p5 = calculate_gravity_term(
        dx1, dy1, dz2, vector_mag_p5, sign_factor=1.0)
    gravity_term_p6 = calculate_gravity_term(
        dx2, dy1, dz2, vector_mag_p6, sign_factor=-1.0)
    gravity_term_p7 = calculate_gravity_term(
        dx1, dy2, dz2, vector_mag_p7, sign_factor=-1.0)
    gravity_term_p8 = calculate_gravity_term(
        dx2, dy2, dz2, vector_mag_p8, sign_factor=1.0)

    delta_gravity = (
        delta_density * gravity_constant * (
            gravity_term_p1
            + gravity_term_p2
            + gravity_term_p3
            + gravity_term_p4
            + gravity_term_p5
            + gravity_term_p6
            + gravity_term_p7
            + gravity_term_p8
        )
    )

    return delta_gravity
end

""" Calculate the gravity term.

# Arguments
- `dx`: x-coordinate of the point relative to the observation point (m)
- `dy`: y-coordinate of the point relative to the observation point (m)
- `dz`: z-coordinate of the point relative to the observation point (m)
- `vector_mag_pt`: Magnitude of the vector from the observation point to the point (m)
- `sign_factor`: Sign factor

# Returns
- `Float64`: Gravity term for a given point
"""
function calculate_gravity_term(
        dx::Float64,
        dy::Float64,
        dz::Float64,
        vector_mag_pt::Float64;
        sign_factor::Float64
)::Float64
    sign_opp = sign_factor * -1.0
    return (
        sign_factor * dz * atan((dx * dy) / (dz * vector_mag_pt))
        + sign_opp * dx * log(vector_mag_pt + dy)
        + sign_opp * dy * log(vector_mag_pt + dx)
    )
end

""" Calculate the magnitude of a vector.
"""
function calculate_vector_magnitude(dx::Float64, dy::Float64, dz::Float64)::Float64
    return sqrt(dx^2 + dy^2 + dz^2)
end

end # module
