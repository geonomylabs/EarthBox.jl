module SolveForVelocity

"""
    solve_for_velocity(
        f1_factor::Float64, f2_factor::Float64, velocity_initial_guess::Float64
    )::Tuple{Float64, Float64}

Solve for conservative orthogonal velocity for top and bottom.

# Arguments
- `f1_factor::Float64`:
  - Factor in the equation for conservative orthogonal velocity equal to:

    2.0*velocity_extension*ysize/xsize (eq. 1)                                     
  
    where velocity_extension is the absolute value of the half extension velocity 
    along the side boundary of the extension model where the magnitude of the 
    extension velocity equal along both side boundaries.

- `f2_factor::Float64`:
  - Factor in the equation for conservative orthogonal velocity equal to:

    (constant1 - 1.0)/constant1 (eq. 2)
    
    where constant1 is the average of the sticky air thickness along the side 
    boundaries of the extension model divided by the ysize of the extension model:
    
    (sticky_thickness_left + sticky_thickness_right)/2.0/ysize (eq. 3)
    
    Note that (eq. 2) simplifies to the following equation  if sticky air 
    thickness is uniform:

    1 - ysize/sticky_air_thickness (eq. 4)

- `velocity_initial_guess::Float64`: 
    - Initial guess for conservative orthogonal velocity.

# Returns
- Tuple containing:
  - `velocity_y_upper::Float64`:
    - Conservative vertical velocity (m/s) along upper boundary of extension model
  - `velocity_y_lower::Float64`:
    - Conservative vertical velocity (m/s) along lower boundary of extension model

"""
function solve_for_velocity(
    f1_factor::Float64,
    f2_factor::Float64,
    velocity_initial_guess::Float64
)::Tuple{Float64, Float64}
    tolerance = 1e-18
    maxiter = 100

    velocity_y_upper_old = velocity_initial_guess
    relative_difference = 1e32
    icount = 0
    while relative_difference > tolerance && icount < maxiter
        if f2_factor == 0.0
            println("f2_factor is zero")
        end
        velocity_y_upper_new = (velocity_y_upper_old - f1_factor)/f2_factor
        if velocity_y_upper_new == 0.0
            println("velocity_y_upper_new is zero")
            println("velocity_y_upper_old ", velocity_y_upper_old)
            println("f1_factor ", f1_factor)
            println("f2_factor ", f2_factor)
        end
        if velocity_y_upper_old == 0.0
            println("velocity_y_upper_old is zero")
            println("velocity_initial_guess ", velocity_initial_guess)
        end

        relative_difference = abs(
            (velocity_y_upper_new - velocity_y_upper_old) / velocity_y_upper_old
        )
        velocity_y_upper_old = velocity_y_upper_new
        icount += 1
    end
    velocity_y_upper = velocity_y_upper_old
    velocity_y_lower = velocity_y_upper*f2_factor
    return (velocity_y_upper, velocity_y_lower)
end

end # module