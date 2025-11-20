module BilinearNumerator

import EarthBox.ModelDataContainer: ModelData
import ..WeightFuncs

""" Calculate marker weight sums for the numerator of the average for heat equation.

This function calculates the sum of weighted marker values for grid nodes, which is used in the numerator
of the bilinear interpolation equation. The marker weight sums are stored in transport arrays with suffix "1".

This function performs this action only for grid scalars associated with the heat equation.

# Arguments
- `model::ModelData`: The model data container containing all necessary arrays and parameters

# Updated Arrays
- `tk1.array::Matrix{Float64}` (ynum, xnum): 
    - temperature (Kelvins) at basic nodes
- `kt1.array::Matrix{Float64}` (ynum, xnum): 
    - conductivity at basic nodes (W/m/K)
- `rhocp1.array::Matrix{Float64}` (ynum, xnum): 
    - rho*Cp on basic nodes
- `hr1.array::Matrix{Float64}` (ynum, xnum): 
    - rad. heat prod. (W/m^3) on basic nodes
- `ha1.array::Matrix{Float64}` (ynum, xnum): 
    - Adiabatic heat production term (expansivity x temperature) on basic grid

# Notes
- Includes serpentinization heat production if the option is active
"""
function calc_marker_weight_sums_for_basic_grid_heat!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    # Number of markers
    marknum = model.markers.parameters.distribution.marknum.value
    
    # Get arrays from model
    #xstp_b = model.grids.arrays.basic.xstp_b.array
    #ystp_b = model.grids.arrays.basic.ystp_b.array
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    
    # Get marker weights
    marker_wtforULnode = model.interpolation.arrays.marker_weights.marker_wtforULnode.array
    marker_wtforLLnode = model.interpolation.arrays.marker_weights.marker_wtforLLnode.array
    marker_wtforURnode = model.interpolation.arrays.marker_weights.marker_wtforURnode.array
    marker_wtforLRnode = model.interpolation.arrays.marker_weights.marker_wtforLRnode.array
    
    # Get marker thermal properties
    marker_TK = model.markers.arrays.thermal.marker_TK.array
    marker_rhocp = model.markers.arrays.thermal.marker_rhocp.array
    marker_kt = model.markers.arrays.thermal.marker_kt.array
    marker_ha = model.markers.arrays.thermal.marker_ha.array
    mat_hr = model.materials.arrays.mat_hr.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    
    # Get output arrays
    tk1 = model.heat_equation.arrays.temperature.tk1.array
    kt1 = model.heat_equation.arrays.thermal_conductivity.kt1.array
    rhocp1 = model.heat_equation.arrays.rhocp.rhocp1.array
    hr1 = model.heat_equation.arrays.radiogenic_production.hr1.array
    ha1 = model.heat_equation.arrays.adiabatic_production.ha1.array
    
    # Get serpentinization properties if needed
    (
        marker_serpentinization_heat_production
    ) = model.markers.arrays.material.marker_serpentinization_heat_production.array
    (
        iuse_serpentinization
    ) = model.materials.parameters.serpentinization.iuse_serpentinization.value
    
    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                # Get indices and weights
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
                
                marker_weight_upr_left = marker_wtforULnode[imarker]
                marker_weight_lower_left = marker_wtforLLnode[imarker]
                marker_weight_upr_right = marker_wtforURnode[imarker]
                marker_weight_lower_right = marker_wtforLRnode[imarker]
                # Get material id
                mitype = marker_matid[imarker]
                mT = marker_TK[imarker]
                mkt = marker_kt[imarker]
                mrhocp = marker_rhocp[imarker]
                radiogenic_heat_production = mat_hr[mitype]
                mha = marker_ha[imarker]
            end

            #@inbounds begin
            #    # Calculate volume factor
            #    vol_fac = WeightFuncs.cell_volume_factor(xstp_b[ix_upr_left], ystp_b[iy_upr_left])
            #end
            
            # Update temperature
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                tk1, mT
            )
            # Update conductivity
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                kt1, mkt
            )
            # Update rho*cp
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                rhocp1, mrhocp
            )
            # Update heat production
            @inbounds begin
                if iuse_serpentinization == 1
                    total_heat_production = radiogenic_heat_production + Float64(marker_serpentinization_heat_production[imarker])
                else
                    total_heat_production = radiogenic_heat_production
                end
            end
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                hr1, total_heat_production
            )
            # Update adiabatic heat production
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                ha1, mha
            )
        end
    end
    
    return nothing
end

""" Check shear modulus for zero value.

# Arguments
- `shear_modulus::Float64`: The shear modulus value to check
- `material_id::Int`: The material ID associated with the shear modulus

# Throws
- `ErrorException`: If the shear modulus is zero

# Notes
- This is a helper function used to validate shear modulus values before calculations
"""
@inline function check_shear_modulus(shear_modulus::Float64, material_id::Int16)::Nothing
    if shear_modulus == 0.0
        error("Shear modulus for marker material is zero. matid = ", material_id)
    end
    return nothing
end

""" Calculate marker weight sums for the numerator of the average for Stokes equations.

This function calculates the sum of weighted marker values for grid nodes, which is used in the numerator
of the bilinear interpolation equation. This function performs this action only for grid scalars associated
with the Stokes-continuity equations.

# Arguments
- `model::ModelData`: The model data container containing all necessary arrays and parameters

# Updated Arrays
- `rho1.array::Matrix{Float64}` (ynum, xnum): 
    - Density (kg/m^3) interpolated from markers to basic grid
- `etan1.array::Matrix{Float64}` (ynum-1, xnum-1): 
    - Normal viscosity (Pa.s) interpolated from markers to pressure grid
- `etas1.array::Matrix{Float64}` (ynum, xnum): 
    - Shear viscosity (Pa.s) interpolated from markers to basic grid
- `eta_flow.array::Matrix{Float64}` (ynum, xnum): 
    - Flow viscosity interpolated from markers to basic grid
- `mus1.array::Matrix{Float64}` (ynum, xnum): 
    - Elastic modulus for shear stress (Pa) interpolated from markers to basic grid
- `mun1.array::Matrix{Float64}` (ynum-1, xnum-1): 
    - Elastic modulus for normal stress (Pa) interpolated from markers to pressure grid
- `sxx1.array::Matrix{Float64}` (ynum-1, xnum-1): 
    - Normal stress (Pa) interpolated from markers to pressure grid
- `sxy1.array::Matrix{Float64}` (ynum, xnum): 
    - Shear stress (Pa) interpolated from markers to basic grid
- `plastics.array::Matrix{Float64}` (ynum, xnum): 
    - Plastic deformation flag interpolated from markers to basic grid
- `plasticn.array::Matrix{Float64}` (ynum-1, xnum-1): 
    - Plastic deformation flag interpolated from markers to pressure grid
- `cohesion_grid.array::Matrix{Float64}` (ynum, xnum): 
    - Cohesion (Pa) interpolated from markers to basic grid
- `fric_degrees_grid.array::Matrix{Float64}` (ynum, xnum): 
    - Friction angle (degrees) interpolated from markers to basic grid
- `dilatation_grid.array::Matrix{Float64}` (ynum, xnum): 
    - Dilatation angle interpolated from markers to basic grid

# Notes
- Harmonic mean is used for shear modulus calculations
- Arithmetic mean is used for all other properties
"""
function calc_marker_weight_sums_for_basic_and_pressure_grids_stokes!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    # Number of markers
    marknum = model.markers.parameters.distribution.marknum.value
    
    # Get arrays from model
    #xstp_b = model.grids.arrays.basic.xstp_b.array
    #ystp_b = model.grids.arrays.basic.ystp_b.array
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array
    
    # Get marker weights
    marker_wtforULnode = model.interpolation.arrays.marker_weights.marker_wtforULnode.array
    marker_wtforLLnode = model.interpolation.arrays.marker_weights.marker_wtforLLnode.array
    marker_wtforURnode = model.interpolation.arrays.marker_weights.marker_wtforURnode.array
    marker_wtforLRnode = model.interpolation.arrays.marker_weights.marker_wtforLRnode.array
    marker_wtforCnode = model.interpolation.arrays.marker_weights.marker_wtforCnode.array
    
    # Get marker properties
    marker_dilatation_angle = model.markers.arrays.rheology.marker_dilatation_angle.array
    marker_rho = model.markers.arrays.material.marker_rho.array
    marker_cohesion = model.markers.arrays.rheology.marker_cohesion.array
    marker_fric = model.markers.arrays.rheology.marker_fric.array
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array
    marker_sxx = model.markers.arrays.stress.marker_sxx.array
    
    marker_sxy = model.markers.arrays.stress.marker_sxy.array

    marker_pfailure = model.markers.arrays.rheology.marker_pfailure.array
    mat_mu = model.materials.arrays.mat_mu.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    # Get output arrays
    rho1 = model.stokes_continuity.arrays.density.rho1.array
    etan1 = model.stokes_continuity.arrays.viscosity.etan1.array
    etas1 = model.stokes_continuity.arrays.viscosity.etas1.array
    eta_flow = model.stokes_continuity.arrays.viscosity.eta_flow.array
    mus1 = model.stokes_continuity.arrays.shear_modulus.mus1.array
    mun1 = model.stokes_continuity.arrays.shear_modulus.mun1.array
    sxx1 = model.stokes_continuity.arrays.stress.sxx1.array
    sxy1 = model.stokes_continuity.arrays.stress.sxy1.array
    plastics = model.stokes_continuity.arrays.plastic_def.plastics.array
    plasticn = model.stokes_continuity.arrays.plastic_def.plasticn.array
    cohesion_grid = model.stokes_continuity.arrays.plastic_def.cohesion_grid.array
    fric_degrees_grid = model.stokes_continuity.arrays.plastic_def.fric_degrees_grid.array
    dilatation_grid = model.stokes_continuity.arrays.plastic_def.dilatation_grid.array
    
    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                # Get indices and weights
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
                dx_upr_left = marker_dx[imarker]
                dy_upr_left = marker_dy[imarker]
                
                marker_weight_upr_left = marker_wtforULnode[imarker]
                marker_weight_lower_left = marker_wtforLLnode[imarker]
                marker_weight_upr_right = marker_wtforURnode[imarker]
                marker_weight_lower_right = marker_wtforLRnode[imarker]
                marker_weight_central = marker_wtforCnode[imarker]
                # Get material id
                mitype = marker_matid[imarker]
                
                mrho = marker_rho[imarker]
                mcohesion = marker_cohesion[imarker]
                mfric = asin(marker_fric[imarker]) * 180.0 / π
                meta_flow = marker_eta_flow[imarker]
                meta = marker_eta[imarker]

                mshear = mat_mu[mitype]
                check_shear_modulus(mshear, mitype)
                mmu_inv = 1/mshear
                
                msxy = marker_sxy[imarker]
                plastyn = Float64(marker_pfailure[imarker])
                mdilatation = Float64(marker_dilatation_angle[imarker])
                msxx = marker_sxx[imarker]
            end
            
            # Calculate volume factor
            #vol_fac = WeightFuncs.cell_volume_factor(xstp_b[ix_upr_left], ystp_b[iy_upr_left])
            
            # Update density
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                rho1, mrho
            )
            # Update cohesion
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                cohesion_grid, mcohesion
            )
            # Update friction angle
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                fric_degrees_grid, mfric
            )
            # Update flow viscosity
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                eta_flow, meta_flow
            )
            # Update viscosity
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                etas1, meta
            )
            # Update normal viscosity
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                etan1, meta
            )
            # Update shear modulus
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                mus1, mmu_inv
            )
            # Update normal shear modulus
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                mun1, mmu_inv
            )
            # Update shear stress
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                sxy1, msxy
            )
            # Update plastic deformation flag
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                plastics, plastyn
            )
            # Update plastic deformation flag for pressure nodes
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                plasticn, plastyn
            )
            # Update dilatation angle
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                dilatation_grid, mdilatation
            )
            # Update normal stress
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                sxx1, msxx
            )
        end
    end
    return nothing
end

""" Calculate marker weight sums for the numerator of the average for vy grid in Stokes equations.

This function calculates the sum of weighted marker values for vy grid nodes, which is used in the numerator
of the bilinear interpolation equation. This function performs this action only for density on the vy grid.

# Arguments
- `model::ModelData`: 
    - The model data container containing all necessary arrays and parameters

# Updated Arrays
- `rho1_vy.array::Matrix{Float64}` (ynum, xnum): 
    - Density (kg/m^3) interpolated from markers to vy grid

"""
function calc_marker_weight_sums_for_vy_grid_stokes!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    # Number of markers
    marknum = model.markers.parameters.distribution.marknum.value
    
    # Get arrays from model
    #xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    #ystp_b = model.grids.arrays.basic.ystp_b.array
    marker_xn_vy = model.markers.arrays.grid_marker_relationship.marker_xn_vy.array
    marker_yn_vy = model.markers.arrays.grid_marker_relationship.marker_yn_vy.array
    
    # Get marker weights
    marker_wtforULnodeVy = model.interpolation.arrays.marker_weights.marker_wtforULnodeVy.array
    marker_wtforLLnodeVy = model.interpolation.arrays.marker_weights.marker_wtforLLnodeVy.array
    marker_wtforURnodeVy = model.interpolation.arrays.marker_weights.marker_wtforURnodeVy.array
    marker_wtforLRnodeVy = model.interpolation.arrays.marker_weights.marker_wtforLRnodeVy.array
    
    # Get output arrays
    rho1_vy = model.stokes_continuity.arrays.density.rho1_vy.array
    marker_rho = model.markers.arrays.material.marker_rho.array
    
    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                # Get indices and weights
                ix_upr_left = marker_xn_vy[imarker]
                iy_upr_left = marker_yn_vy[imarker]
                
                marker_weight_upr_left = marker_wtforULnodeVy[imarker]
                marker_weight_lower_left = marker_wtforLLnodeVy[imarker]
                marker_weight_upr_right = marker_wtforURnodeVy[imarker]
                marker_weight_lower_right = marker_wtforLRnodeVy[imarker]

                mrho = marker_rho[imarker]
            end
            # Calculate volume factor
            #vol_fac = WeightFuncs.cell_volume_factor(xstp_vy[ix_upr_left], ystp_b[iy_upr_left])
            WeightFuncs.update_weight_inclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                rho1_vy, mrho
            )
        end
    end
    
    return nothing
end

""" Calculate marker weight sums for viscosity.

This function calculates the sum of weighted marker values for viscosity on both basic and pressure grids.

# Arguments
- `model::ModelData`: The model data container containing all necessary arrays and parameters

# Updated Arrays
- `etan1.array::Matrix{Float64}` (ynum-1, xnum-1):
    - Normal viscosity (Pa.s) interpolated from markers to pressure grid
- `etas1.array::Matrix{Float64}` (ynum, xnum):
    - Shear viscosity (Pa.s) interpolated from markers to basic grid

"""
function calc_marker_weight_sums_viscosity!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    # Number of markers
    marknum = model.markers.parameters.distribution.marknum.value
    
    # Get arrays from model
    #xstp_b = model.grids.arrays.basic.xstp_b.array
    #ystp_b = model.grids.arrays.basic.ystp_b.array
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array
    
    # Get marker weights
    marker_wtforULnode = model.interpolation.arrays.marker_weights.marker_wtforULnode.array
    marker_wtforLLnode = model.interpolation.arrays.marker_weights.marker_wtforLLnode.array
    marker_wtforURnode = model.interpolation.arrays.marker_weights.marker_wtforURnode.array
    marker_wtforLRnode = model.interpolation.arrays.marker_weights.marker_wtforLRnode.array
    marker_wtforCnode = model.interpolation.arrays.marker_weights.marker_wtforCnode.array
    
    # Get marker properties
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    
    # Get output arrays
    etan1 = model.stokes_continuity.arrays.viscosity.etan1.array
    etas1 = model.stokes_continuity.arrays.viscosity.etas1.array
    
    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                # Get indices and weights
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
                dx_upr_left = marker_dx[imarker]
                dy_upr_left = marker_dy[imarker]
                
                marker_weight_upr_left = marker_wtforULnode[imarker]
                marker_weight_lower_left = marker_wtforLLnode[imarker]
                marker_weight_upr_right = marker_wtforURnode[imarker]
                marker_weight_lower_right = marker_wtforLRnode[imarker]
                marker_weight_central = marker_wtforCnode[imarker]

                meta = marker_eta[imarker]
            end
            # Calculate volume factor
            #vol_fac = WeightFuncs.cell_volume_factor(xstp_b[ix_upr_left], ystp_b[iy_upr_left])
            # Update shear viscosity
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                etas1, meta
            )
            
            # Update normal viscosity
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                etan1, meta
            )
        end
    end
    
    return nothing
end

end # module 