module BilinearAverage

import EarthBox.ModelDataContainer.ModelData

"""
    calc_marker_average_at_basic_nodes_heat!(model::ModelData)::Nothing

Calculate 1st order bilinear average marker values at nodes for heat equation.
"""
function calc_marker_average_at_basic_nodes_heat!(model::ModelData)::Nothing
    # Sum of marker weights at basic nodes using an inclusive search radius
    sum_of_marker_weights_at_nodes_basic_inclusive = model.interpolation.arrays.grid_weights.wtnodes.array
    
    # Temperature on basic grid
    calc_bilinear_average_arithmetic!(
        model.heat_equation.arrays.temperature.tk1.array,
        model.heat_equation.arrays.temperature.tk0.array,
        sum_of_marker_weights_at_nodes_basic_inclusive
    )
    
    # Thermal conductivity on basic grid
    calc_bilinear_average_arithmetic!(
        model.heat_equation.arrays.thermal_conductivity.kt1.array,
        model.heat_equation.arrays.thermal_conductivity.kt0.array,
        sum_of_marker_weights_at_nodes_basic_inclusive
    )
    
    # Density times heat capacity on basic grid
    calc_bilinear_average_arithmetic!(
        model.heat_equation.arrays.rhocp.rhocp1.array,
        model.heat_equation.arrays.rhocp.rhocp0.array,
        sum_of_marker_weights_at_nodes_basic_inclusive
    )
    
    # Radiogenic heat production on basic grid
    calc_bilinear_average_arithmetic!(
        model.heat_equation.arrays.radiogenic_production.hr1.array,
        model.heat_equation.arrays.radiogenic_production.hr0.array,
        sum_of_marker_weights_at_nodes_basic_inclusive
    )
    
    # Adiabatic heat production term on basic grid
    calc_bilinear_average_arithmetic!(
        model.heat_equation.arrays.adiabatic_production.ha1.array,
        model.heat_equation.arrays.adiabatic_production.ha0.array,
        sum_of_marker_weights_at_nodes_basic_inclusive
    )
    
    return nothing
end

"""
    calc_marker_average_at_staggered_nodes_stokes!(model::ModelData)::Nothing

Calculate 1st order bilinear average marker values at nodes for Stokes equations.
"""
function calc_marker_average_at_staggered_nodes_stokes!(model::ModelData)::Nothing
    # Sum of marker weights at basic nodes using an inclusive search radius
    sum_of_marker_weights_at_nodes_basic_inclusive = model.interpolation.arrays.grid_weights.wtnodes.array
    
    # Density on basic grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.density.rho1.array,
        model.stokes_continuity.arrays.density.rho0.array,
        sum_of_marker_weights_at_nodes_basic_inclusive
    )
    
    # Sum of marker weights at vy nodes using an inclusive search radius
    sum_of_marker_weights_at_nodes_vy_inclusive = model.interpolation.arrays.grid_weights.wtnodes_vy.array
    
    # Density on vy grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.density.rho1_vy.array,
        model.stokes_continuity.arrays.density.rho0_vy.array,
        sum_of_marker_weights_at_nodes_vy_inclusive
    )
    
    # Sum of marker weights at basic nodes using an exclusive search radius
    sum_of_marker_weights_at_nodes_basic_exclusive = model.interpolation.arrays.grid_weights.wtetas.array
    
    # Shear modulus on basic grid
    calc_bilinear_average_harmonic!(
        model.stokes_continuity.arrays.shear_modulus.mus1.array,
        model.stokes_continuity.arrays.shear_modulus.mus0.array,
        sum_of_marker_weights_at_nodes_basic_exclusive
    )
    
    # Shear viscosity on basic grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.viscosity.etas1.array,
        model.stokes_continuity.arrays.viscosity.etas0.array,
        sum_of_marker_weights_at_nodes_basic_exclusive
    )
    
    check_2d_array_for_zeros!("etas1", model.stokes_continuity.arrays.viscosity.etas1.array)
    
    # Shear stress on basic grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.stress.sxy1.array,
        model.stokes_continuity.arrays.stress.sxy0.array,
        sum_of_marker_weights_at_nodes_basic_exclusive
    )

    # Shear plastic deformation flag on basic grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.plastic_def.plastics.array,
        model.stokes_continuity.arrays.plastic_def.plastics0.array,
        sum_of_marker_weights_at_nodes_basic_exclusive
    )
    
    # Cohesion on basic grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.plastic_def.cohesion_grid.array,
        model.stokes_continuity.arrays.plastic_def.cohesion_grid0.array,
        sum_of_marker_weights_at_nodes_basic_exclusive
    )
    
    # Friction angle on basic grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.plastic_def.fric_degrees_grid.array,
        model.stokes_continuity.arrays.plastic_def.fric_degrees_grid0.array,
        sum_of_marker_weights_at_nodes_basic_exclusive
    )
    
    # Flow viscosity on basic grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.viscosity.eta_flow.array,
        model.stokes_continuity.arrays.viscosity.eta_flow0.array,
        sum_of_marker_weights_at_nodes_basic_exclusive
    )

    # Sum of marker weights at pressure nodes using an inclusive search radius
    sum_of_marker_weights_at_nodes_pressure_inclusive = model.interpolation.arrays.grid_weights.wtetan.array
    
    # Shear modulus on pressure grid
    calc_bilinear_average_harmonic!(
        model.stokes_continuity.arrays.shear_modulus.mun1.array,
        model.stokes_continuity.arrays.shear_modulus.mun0.array,
        sum_of_marker_weights_at_nodes_pressure_inclusive
    )
    
    # Normal viscosity on pressure grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.viscosity.etan1.array,
        model.stokes_continuity.arrays.viscosity.etan0.array,
        sum_of_marker_weights_at_nodes_pressure_inclusive
    )
    
    # Normal stress on pressure grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.stress.sxx1.array,
        model.stokes_continuity.arrays.stress.sxx0.array,
        sum_of_marker_weights_at_nodes_pressure_inclusive
    )
    
    # Normal plastic deformation flag on pressure grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.plastic_def.plasticn.array,
        model.stokes_continuity.arrays.plastic_def.plasticn0.array,
        sum_of_marker_weights_at_nodes_pressure_inclusive
    )
    
    # Dilatation on pressure grid
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.plastic_def.dilatation_grid.array,
        model.stokes_continuity.arrays.plastic_def.dilatation_grid0.array,
        sum_of_marker_weights_at_nodes_pressure_inclusive
    )
    
    return nothing
end

"""
    check_2d_array_for_zeros(array_name::String, array::Matrix{Float64})::Nothing

Check for zero values in a 2D array and print warnings.
"""
function check_2d_array_for_zeros!(array_name::String, array::Matrix{Float64})::Nothing
    ynum, xnum = size(array)
    for j in 1:xnum
        for i in 1:ynum
            if array[i, j] == 0
                print_zero_warning(array_name, i, j)
            end
        end
    end
    return nothing
end

"""
    print_zero_warning(array_name::String, i::Int, j::Int)::Nothing

Print warning message if zero value is found in array.
"""
function print_zero_warning(array_name::String, i::Int, j::Int)::Nothing
    println("For 2D array ", array_name, ", a zero value was found at index i = ", i, 
            " j = ", j, ". Consider increasing marker resolution.")
    return nothing
end

"""
    calc_marker_average_at_nodes_viscosity!(model::ModelData)::Nothing

Calculate average marker values at nodes for viscosity.
"""
function calc_marker_average_at_nodes_viscosity!(model::ModelData)::Nothing
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.viscosity.etas1.array,
        model.stokes_continuity.arrays.viscosity.etas0.array,
        model.interpolation.arrays.grid_weights.wtetas.array
    )
    
    calc_bilinear_average_arithmetic!(
        model.stokes_continuity.arrays.viscosity.etan1.array,
        model.stokes_continuity.arrays.viscosity.etan0.array,
        model.interpolation.arrays.grid_weights.wtetan.array
    )
    
    return nothing
end

"""
    calc_marker_average_subgrid_stress_change!(model::ModelData)::Nothing

Finalize grid stress change by dividing by integrated marker weights.
"""
function calc_marker_average_subgrid_stress_change!(model::ModelData)::Nothing
    grid_weights = model.interpolation.arrays.grid_weights
    stress_change = model.stokes_continuity.arrays.stress_change
    
    calc_bilinear_average_arithmetic!(
        stress_change.dsxyn.array,
        copy(stress_change.dsxyn.array),
        grid_weights.wtetas.array
    )
    
    calc_bilinear_average_arithmetic!(
        stress_change.dsxxn.array,
        copy(stress_change.dsxxn.array),
        grid_weights.wtetan.array
    )
    
    return nothing
end

"""
    calc_marker_average_subgrid_heat_diffusion!(model::ModelData)::Nothing

Finalize grid stress change by dividing by integrated marker weights.
"""
function calc_marker_average_subgrid_heat_diffusion!(model::ModelData)::Nothing
    calc_bilinear_average_arithmetic!(
        model.heat_equation.arrays.temperature.dtkn.array,
        copy(model.heat_equation.arrays.temperature.dtkn.array),
        model.interpolation.arrays.grid_weights.wtnodes.array
    )
    
    return nothing
end

"""
    calc_bilinear_average_arithmetic!(
        bilinear_avg::Matrix{Float64}, bilinear_avg_old::Matrix{Float64}, 
        sum_of_marker_weights_at_nodes::Matrix{Float64}
    )::Nothing

Calculate bilinear average using arithmetic mean at grid nodes.
"""
function calc_bilinear_average_arithmetic!(
    bilinear_avg::Matrix{Float64}, 
    bilinear_avg_old::Matrix{Float64}, 
    sum_of_marker_weights_at_nodes::Matrix{Float64}
)::Nothing
    ynum, xnum = size(bilinear_avg)
    Threads.@threads for j in 1:xnum
        for i in 1:ynum
            @inbounds begin
                if sum_of_marker_weights_at_nodes[i, j] != 0
                    bilinear_avg[i, j] = bilinear_avg[i, j] / sum_of_marker_weights_at_nodes[i, j]
                else
                    bilinear_avg[i, j] = bilinear_avg_old[i, j]
                end
            end
        end
    end
    return nothing
end

"""
    calc_bilinear_average_harmonic!(
        bilinear_avg::Matrix{Float64}, bilinear_avg_old::Matrix{Float64}, 
        sum_of_marker_weights_at_nodes::Matrix{Float64}
    )::Nothing

Calculate bilinear average using harmonic mean at grid nodes.
"""
function calc_bilinear_average_harmonic!(
    bilinear_avg::Matrix{Float64}, 
    bilinear_avg_old::Matrix{Float64}, 
    sum_of_marker_weights_at_nodes::Matrix{Float64}
)::Nothing
    ynum, xnum = size(bilinear_avg)
    Threads.@threads for j in 1:xnum
        for i in 1:ynum
            @inbounds begin
                if sum_of_marker_weights_at_nodes[i, j] != 0
                    bilinear_avg[i, j] = 1.0 / (bilinear_avg[i, j] / sum_of_marker_weights_at_nodes[i, j])
                else
                    bilinear_avg[i, j] = bilinear_avg_old[i, j]
                end
            end
        end
    end
    return nothing
end

end # module 