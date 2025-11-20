module HeatResiduals
import EarthBox.ModelDataContainer: ModelData
import EarthBox.Domain: basic_node_on_boundary

""" Compute residuals for heat equation solution.

Temperature equation: RHO*CP*DT/Dt = -dqx/dx - dqy/dy + Ht

# Updated Arrays
- `model.heat_equation.arrays.residuals.rest.array`: Temperature residuals
"""
function compute_residuals_temp!(model::ModelData)
    compute_residual_on_boundary!(model)
    compute_residual_time_derivative_term!(model)
    compute_residual_heat_flux_x_derivative_term!(model)
    compute_residual_heat_flux_y_derivative_term!(model)
end

""" Compute residuals for heat equation solution on boundaries.

Temperature equation: RHO*CP*DT/Dt = -dqx/dx - dqy/dy + Ht

# Updated Arrays
- `model.heat_equation.arrays.residuals.rest.array`: Temperature residuals
"""
function compute_residual_on_boundary!(model::ModelData)
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    rest = model.heat_equation.arrays.residuals.rest.array
    
    for j in 1:xnum
        for i in 1:ynum
            if basic_node_on_boundary(i, j, xnum, ynum)
                rest[i, j] = 0.0
            end
        end
    end
end

""" Compute residual for term with time derivative Ht-DT/dt.

Temperature equation: RHO*CP*DT/Dt = -dqx/dx - dqy/dy + Ht

# Updated Arrays
- `model.heat_equation.arrays.residuals.rest.array`: Temperature residuals
"""
function compute_residual_time_derivative_term!(model::ModelData)
    timestep = model.timestep.parameters.thermal_loop.timestep_heat.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    RT = model.heat_equation.arrays.rhs.RT1.array
    rhocp = model.heat_equation.arrays.rhocp.rhocp1.array
    tknew = model.heat_equation.arrays.temperature.tk2.array
    tk = model.heat_equation.arrays.temperature.tk0.array
    rest = model.heat_equation.arrays.residuals.rest.array
    
    for j in 1:xnum
        for i in 1:ynum
            if !basic_node_on_boundary(i, j, xnum, ynum)
                rest[i, j] = RT[i, j] + time_derivative_term(
                    (i, j), rhocp, tk, tknew, timestep)
            end
        end
    end
end

""" Time derivative term for heat residual: -DT/dt.
"""
@inline function time_derivative_term(
    ij_indices::Tuple{Int,Int}, 
    rhocp::Matrix{Float64},
    tk::Matrix{Float64}, 
    tknew::Matrix{Float64},
    timestep::Float64
)::Float64
    i, j = ij_indices
    return -rhocp[i, j] * (tknew[i, j] - tk[i, j]) / timestep
end

""" Compute residuals for term with heat flux x derivative -dqx/dx.

Temperature equation: RHO*CP*DT/Dt = -dqx/dx - dqy/dy + Ht

# Updated Arrays
- `model.heat_equation.arrays.residuals.rest.array`: Temperature residuals
"""
function compute_residual_heat_flux_x_derivative_term!(model::ModelData)
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xstp = model.grids.arrays.basic.xstp_b.array
    xstpc = model.grids.arrays.pressure.xstp_pr.array
    tknew = model.heat_equation.arrays.temperature.tk2.array
    kt = model.heat_equation.arrays.thermal_conductivity.kt1.array
    rest = model.heat_equation.arrays.residuals.rest.array
    
    for j in 1:xnum
        for i in 1:ynum
            if !basic_node_on_boundary(i, j, xnum, ynum)
                rest[i, j] += heat_flux_x_derivative_term(
                    (i, j), xstp, xstpc, kt, tknew)
            end
        end
    end
end

""" Heat flux x derivative term for heat residual.
"""
@inline function heat_flux_x_derivative_term(
    ij_indices::Tuple{Int,Int},
    xstp::Vector{Float64}, 
    xstpc::Vector{Float64},
    kt::Matrix{Float64}, 
    tknew::Matrix{Float64}
)::Float64
    i, j = ij_indices
    return ((kt[i, j] + kt[i, j+1]) * (tknew[i, j+1] - tknew[i, j]) / xstp[j] -
            (kt[i, j-1] + kt[i, j]) * (tknew[i, j] - tknew[i, j-1]) / xstp[j-1]) /
           2.0 / xstpc[j-1]
end

""" Compute residuals for term with heat flux y derivative -dqy/dy.

Temperature equation: RHO*CP*DT/Dt = -dqx/dx - dqy/dy + Ht

# Updated Arrays
- `model.heat_equation.arrays.residuals.rest.array`: Temperature residuals
"""
function compute_residual_heat_flux_y_derivative_term!(model::ModelData)
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    ystp = model.grids.arrays.basic.ystp_b.array
    ystpc = model.grids.arrays.pressure.ystp_pr.array
    tknew = model.heat_equation.arrays.temperature.tk2.array
    kt = model.heat_equation.arrays.thermal_conductivity.kt1.array
    rest = model.heat_equation.arrays.residuals.rest.array
    
    for j in 1:xnum
        for i in 1:ynum
            if !basic_node_on_boundary(i, j, xnum, ynum)
                rest[i, j] += heat_flux_y_derivative_term(
                    (i, j), ystp, ystpc, kt, tknew)
            end
        end
    end
end

""" Heat flux y derivative term for heat residual.
"""
@inline function heat_flux_y_derivative_term(
    ij_indices::Tuple{Int,Int},
    ystp::Vector{Float64}, 
    ystpc::Vector{Float64},
    kt::Matrix{Float64}, 
    tknew::Matrix{Float64}
)::Float64
    i, j = ij_indices
    return ((kt[i, j] + kt[i+1, j]) * (tknew[i+1, j] - tknew[i, j]) / ystp[i] -
            (kt[i-1, j] + kt[i, j]) * (tknew[i, j] - tknew[i-1, j]) / ystp[i-1]) /
           2.0 / ystpc[i-1]
end

end # module
