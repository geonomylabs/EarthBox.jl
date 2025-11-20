module PlasticFailureNodes

import EarthBox.ModelDataContainer: ModelData
import EarthBox: MobileWall
import ..PlasticFailure: CheckYield, YieldStress

"""
    update_nodes_for_plastic_yielding(model, no_yielding_in_mobile_wall, no_yielding_in_plate_extension)

Update basic nodes for plastic yielding.

# Arguments
- `model`: Model data structure
- `no_yielding_in_mobile_wall`: Boolean to prevent yielding in mobile wall
- `no_yielding_in_plate_extension`: Boolean to prevent yielding in plate extension

# Returns
- `etavp_new`: New viscoplastic viscosity on basic grid
- `plastic_strain_rate`: Plastic strain rate
- `plastic_yield_new`: New yielding indicators on basic grid
- `nnode_yield`: Number of basic nodes that have reached yield criteria
- `global_yield_error`: Global yield error
"""
function update_nodes_for_plastic_yielding(
    model::ModelData,
    no_yielding_in_mobile_wall::Bool,
    no_yielding_in_plate_extension::Bool = false
)::Tuple{
    Matrix{Float64}, # etavp_new
    Matrix{Float64}, # plastic_strain_rate
    Matrix{Float64}, # plastic_yield_new
    Int64, # nnode_yield
    Float64 # global_yield_error
}
    # General parameters
    timestep          = model.timestep.parameters.main_time_loop.timestep.value
    timesum           = model.timestep.parameters.main_time_loop.timesum.value
    xnum              = model.grids.parameters.geometry.xnum.value
    ynum              = model.grids.parameters.geometry.ynum.value
    viscosity_minimum = model.materials.parameters.viscosity_limits.viscosity_min.value
    viscosity_maximum = model.materials.parameters.viscosity_limits.viscosity_max.value
    # Fluid pressure model parameters
    iuse_fluid_pressure_for_yield = model.materials.parameters.stress_limits_yield.iuse_fluid_pressure_for_yield.value
    ntimestep                     = model.timestep.parameters.main_time_loop.ntimestep.value
    y_sealevel                    = model.topography.parameters.sealevel.y_sealevel.value
    # Mobile wall parameters
    velocity_internal_x       = model.bcs.parameters.velocity.velocity_internal_x.value
    x_left_mobile_wall        = model.geometry.parameters.mobile_wall.x_left_mobile_wall.value
    x_right_mobile_wall       = model.geometry.parameters.mobile_wall.x_right_mobile_wall.value
    y_top_mobile_wall         = model.geometry.parameters.mobile_wall.y_top_mobile_wall.value
    y_bottom_mobile_wall      = model.geometry.parameters.mobile_wall.y_bottom_mobile_wall.value
    plate_extension_thickness = model.geometry.parameters.mobile_wall.plate_extension_thickness.value
    plate_extension_width     = model.geometry.parameters.mobile_wall.plate_extension_width.value
    # Grid scalar arrays
    cohesion_grid     = model.stokes_continuity.arrays.plastic_def.cohesion_grid.array
    fric_degrees_grid = model.stokes_continuity.arrays.plastic_def.fric_degrees_grid.array
    plastic_yield     = model.stokes_continuity.arrays.plastic_def.plastic_yield.array
    yield_error       = model.stokes_continuity.arrays.plastic_def.yield_error.array
    etas1             = model.stokes_continuity.arrays.viscosity.etas1.array
    eta_flow          = model.stokes_continuity.arrays.viscosity.eta_flow.array
    sxy2              = model.stokes_continuity.arrays.stress.sxy2.array
    sxx2              = model.stokes_continuity.arrays.stress.sxx2.array
    mus1              = model.stokes_continuity.arrays.shear_modulus.mus1.array
    pr1               = model.stokes_continuity.arrays.pressure.pr1.array
    # 1D basic grid arrays
    gridx = model.grids.arrays.basic.gridx_b.array
    gridy = model.grids.arrays.basic.gridy_b.array

    viscosity_weight_factor = 0.0

    etavp_new = copy(eta_flow)
    plastic_yield_new = zeros(ynum, xnum)
    plastic_strain_rate = zeros(ynum, xnum)
    nnode_yield = 0
    global_yield_error = 0.0
    global_yield_error_sum = 0.0
    pi_term = π / 180.0

    for j in 1:xnum
        for i in 1:ynum
            yield_marker               = plastic_yield[i, j]
            viscosity_viscoplastic_old = etas1[i, j]
            viscosity_flow             = eta_flow[i, j]
            shear_modulus              = mus1[i, j]
            cohesion                   = cohesion_grid[i, j]
            friction_coef              = sin(fric_degrees_grid[i, j] * pi_term)

            yield_is_applicable = CheckYield.check_yield_is_applicable(
                ntimestep, cohesion, friction_coef)
            # Ensure that yielding does not occur in nodes located within the
            # mobile wall or plate extension from mobile wall
            if no_yielding_in_mobile_wall
                in_wall = MobileWall.in_mobile_wall_transient(
                    gridx[j], gridy[i],
                    x_left_mobile_wall, x_right_mobile_wall,
                    y_top_mobile_wall, y_bottom_mobile_wall,
                    velocity_internal_x, timesum
                )
                if in_wall
                    yield_is_applicable = false
                end
            end
            if no_yielding_in_plate_extension
                in_plate_extension = MobileWall.in_plate_extension_transient(
                    gridx[j], gridy[i], plate_extension_thickness,
                    plate_extension_width, x_left_mobile_wall,
                    y_bottom_mobile_wall, velocity_internal_x, timesum
                )
                if in_plate_extension
                    yield_is_applicable = false
                end
            end
            if yield_is_applicable
                pressure = interpolate_from_pressure_grid_to_basic_grid(i, j, ynum, xnum, pr1)
                stress_invariant = calc_stress_invariant_on_grid(i, j, ynum, xnum, sxy2, sxx2)
                stress_invariant_elastic = calc_stress_invariant_elastic(
                    stress_invariant, shear_modulus,
                    viscosity_viscoplastic_old, timestep
                )
                stress_yield = YieldStress.calc_yield_stress(
                    cohesion, friction_coef, pressure,
                    iuse_fluid_pressure_for_yield, gridy[i], y_sealevel)
                node_at_yield = check_for_yielding(yield_marker)
                if node_at_yield
                    stress_error = calc_nodal_stress_error_at_yield(stress_invariant, stress_yield)
                    global_yield_error_sum = update_yield_error_sum(
                        global_yield_error_sum, stress_error)
                    yield_error[i, j] = stress_error
                    nnode_yield += 1
                end
                if stress_yield < stress_invariant_elastic
                    viscosity_viscoplastic_new = calc_viscosity_at_yield(
                        timestep, shear_modulus, stress_yield, stress_invariant_elastic)
                    if viscosity_viscoplastic_new < viscosity_flow
                        plastic_yield_new[i, j] = 1
                        viscosity_viscoplastic_new = calc_weighted_viscosity(
                            viscosity_viscoplastic_old, viscosity_viscoplastic_new,
                            viscosity_weight_factor
                        )
                        viscosity_viscoplastic_new = apply_viscosity_limits(
                            viscosity_minimum, viscosity_maximum, viscosity_viscoplastic_new)
                        etavp_new[i, j] = viscosity_viscoplastic_new
                        plastic_strain_rate[i, j] = (
                            0.5*stress_invariant/viscosity_viscoplastic_new
                            - 0.5*stress_invariant/viscosity_flow
                        )
                        if !node_at_yield
                            stress_error = calc_nodal_stress_error_at_yield(
                                stress_invariant, stress_yield)
                            global_yield_error_sum = update_yield_error_sum(
                                global_yield_error_sum, stress_error)
                            yield_error[i, j] = stress_error
                            nnode_yield += 1
                        end
                    else
                        etavp_new[i, j] = viscosity_flow
                        plastic_yield_new[i, j] = 0
                    end
                else
                    etavp_new[i, j] = viscosity_flow
                    plastic_yield_new[i, j] = 0
                end
            else
                etavp_new[i, j] = viscosity_flow
                plastic_yield_new[i, j] = 0
            end
        end
    end
    if nnode_yield > 0
        global_yield_error = calc_global_yield_error(global_yield_error_sum, nnode_yield)
    end
    return (
        etavp_new, plastic_strain_rate, plastic_yield_new,
        nnode_yield, global_yield_error
    )
end

@inline function calc_stress_invariant_on_grid(
    i::Int,
    j::Int,
    ynum::Int,
    xnum::Int,
    sxy2::Matrix{Float64},
    sxx2::Matrix{Float64}
)::Float64
    sxx = interpolate_from_pressure_grid_to_basic_grid(
        i, j, ynum, xnum, sxx2)
    sxy = sxy2[i, j]
    stress_invariant = sqrt(sxy^2.0 + sxx^2.0)
    return stress_invariant
end

@inline function interpolate_from_pressure_grid_to_basic_grid(
    i::Int,
    j::Int,
    ynum::Int,
    xnum::Int,
    pr_array::Matrix{Float64}
)::Float64
    scalar = 0.0
    # Internal basic nodes
    if (1 < i < ynum && 1 < j < xnum)
        scalar = (
            pr_array[i-1, j-1] +
            pr_array[i, j-1] +
            pr_array[i-1, j] +
            pr_array[i, j]
        ) / 4.0
    # Corner nodes
    elseif i == 1 && j == 1
        scalar = pr_array[i, j]
    elseif i == ynum && j == 1
        scalar = pr_array[i-1, j]
    elseif i == ynum && j == xnum
        scalar = pr_array[i-1, j-1]
    elseif i == 1 && j == xnum
        scalar = pr_array[i, j-1]
    # Left boundary
    elseif j == 1
        scalar = (pr_array[i, j] + pr_array[i-1, j]) / 2.0
    # Right boundary
    elseif j == xnum
        scalar = (pr_array[i-1, j-1] + pr_array[i, j-1]) / 2.0
    # Top boundary
    elseif i == 1
        scalar = (pr_array[i, j] + pr_array[i, j-1]) / 2.0
    # Bottom boundary
    elseif i == ynum
        scalar = (pr_array[i-1, j] + pr_array[i-1, j-1]) / 2.0
    end
    return scalar
end

@inline function calc_stress_invariant_elastic(
    stress_invariant::Float64,
    shear_modulus::Float64,
    viscosity_viscoplastic_old::Float64,
    timestep::Float64
)::Float64
    stress_invariant_elastic = (
        stress_invariant
        * (shear_modulus * timestep + viscosity_viscoplastic_old)
        / viscosity_viscoplastic_old
    )
    return stress_invariant_elastic
end

@inline function calc_nodal_stress_error_at_yield(
    stress_invariant::Float64,
    stress_yield::Float64
)::Float64
    stress_error = stress_invariant - stress_yield
    return stress_error
end

@inline function update_yield_error_sum(
    global_yield_error_sum::Float64,
    stress_error::Float64
)::Float64
    global_yield_error_sum = global_yield_error_sum + stress_error^2.0
    return global_yield_error_sum
end

@inline function calc_viscosity_at_yield(
    timestep::Float64,
    shear_modulus::Float64,
    stress_yield::Float64,
    stress_invariant_elastic::Float64
)::Float64
    viscosity = (
        timestep
        * shear_modulus
        * stress_yield / (stress_invariant_elastic - stress_yield)
    )
    return viscosity
end

@inline function calc_weighted_viscosity(
    viscosity_viscoplastic_old::Float64,
    viscosity_viscoplastic_new::Float64,
    viscosity_weight_factor::Float64
)::Float64
    viscosity_viscoplastic_new = (
        viscosity_viscoplastic_new^(1.0 - viscosity_weight_factor)
        * viscosity_viscoplastic_old^viscosity_weight_factor
    )
    return viscosity_viscoplastic_new
end

@inline function check_for_yielding(yield_marker::Float64)::Bool
    return yield_marker > 0
end

@inline function apply_viscosity_limits(
    viscosity_minimum::Float64,
    viscosity_maximum::Float64,
    viscosity::Float64
)::Float64
    viscosity = max(viscosity, viscosity_minimum)
    viscosity = min(viscosity, viscosity_maximum)
    return viscosity
end

@inline function calc_global_yield_error(
    global_yield_error_sum::Float64,
    nnode_yield::Int
)::Float64
    global_yield_error = sqrt(global_yield_error_sum / nnode_yield)
    return global_yield_error
end

end # module 