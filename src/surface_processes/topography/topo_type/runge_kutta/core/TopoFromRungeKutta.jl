module TopoFromRungeKutta

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Advection: MarkersLocation
import EarthBox.Advection: VelocityAndSpin
import EarthBox: MathTools

""" Calculate topography from markers and grid velocity.

Topography markers are advected to new locations using the Runge-Kutta method and new 
elevations are interpolated back to the master topography grid with fixed locations. Note 
that this function does not account for erosion and depositional processes.

# Arguments
- `model::ModelData`: The model data containing topography information
- `runge_kutta_order_max::Int`: Maximum order of Runge-Kutta method to use

# Updated Arrays
- `model.topography.arrays.gridt.array`: Array((7, toponum), Float64)
    The sub-array gridt[2,toponum], which stores the y-coordinate of topography is updated.
"""
function calc_topo_from_runge_kutta!(
    model::ModelData,
    inside_flags::Vector{Int8},
    runge_kutta_order_max::Int
)::Nothing
    timestep = model.timestep.parameters.main_time_loop.timestep.value

    gridt = model.topography.arrays.gridt.array
    topo_gridx = copy(gridt[1, :])
    topo_gridy = copy(gridt[2, :])
    toponum = length(topo_gridx)

    (
        topo_velocity_x,
        topo_velocity_y
    ) = VelocityAndSpin.interpolate_velocity_for_topography_markers(
        model, topo_gridx, topo_gridy, inside_flags, runge_kutta_order_max)

    topo_gridx_advect = zeros(Float64, toponum)
    topo_gridy_advect = zeros(Float64, toponum)

    for i in 1:toponum
        topo_gridx_advect[i] = MarkersLocation.advect_marker(
            topo_gridx[i], topo_velocity_x[i], timestep)
        topo_gridy_advect[i] = MarkersLocation.advect_marker(
            topo_gridy[i], topo_velocity_y[i], timestep)
    end
    topo_gridy_new = zeros(Float64, toponum)

    MathTools.linear_interp_vals!(
        topo_gridx_advect,
        topo_gridy_advect,
        topo_gridx,
        topo_gridy_new
    )

    for i in 1:toponum
        gridt[2, i] = topo_gridy_new[i]
    end

    return nothing
end

""" Print information about topography update for debugging purposes.
"""
function print_info(
    toponum::Int,
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_gridx_advect::Vector{Float64},
    topo_gridy_advect::Vector{Float64},
    topo_gridy_new::Vector{Float64},
    topo_velocity_x::Vector{Float64},
    topo_velocity_y::Vector{Float64},
    gridt::Matrix{Float64}
)::Nothing
    println("     >> Topography Update Information")
    i = 701  # Adjusted for 1-based indexing
    xo = round(topo_gridx[i], digits=2)
    yo = round(topo_gridy[i], digits=2)
    xad = round(topo_gridx_advect[i], digits=2)
    yad = round(topo_gridy_advect[i], digits=2)
    xint = round(topo_gridx[i], digits=2)
    yint = round(topo_gridy_new[i], digits=2)
    con_fac_to_cm_yr = 100.0 * 365.25 * 24.0 * 3600.0
    vxrk = round(topo_velocity_x[i] * con_fac_to_cm_yr, digits=5)
    vyrk = round(topo_velocity_y[i] * con_fac_to_cm_yr, digits=5)
    vx = round(gridt[4, i] * con_fac_to_cm_yr, digits=5)
    vy = round(gridt[5, i] * con_fac_to_cm_yr, digits=5)
    println("     >> i xo yo xad yad xint yint")
    println("     >> TOPO UPDATE: ", i, " ", xo, " ", yo, " ", xad, " ", yad, " ", xint, " ", 
            yint)
    println("     >> i vxrk, vyrk, vx, vy")
    println("     >> TOPO VELOC: ", i, " ", vxrk, " ", vyrk, " ", vx, " ", vy)
    
    return nothing
end

end # module 