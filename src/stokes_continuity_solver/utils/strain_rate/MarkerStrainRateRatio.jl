module MarkerStrainRateRatio

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Interpolation: MarkerGridMapping
import EarthBox.Interpolation: GridToMarker
import EarthBox: GridFuncs

""" Call function to calculate marker-grid strain rate ratio for markers.

The strain rate ratio is used in the plasticity model and composite viscosity 
model to convert velocity-based strain rate invariants interpolated from grid to 
strain-rate invariants based on nodal stress. As described in Gerya (2010), the 
use of nodal stress to calculate strain rate invariants reduces numerical 
diffusion in shear zones.

The strain rate ratio is given by the following equation:

    marker_sr_ratio = eii_stress/eii_velocity                          eq. 1

where eii_velocity is the grid-velocity based second invariant of strain rate 
invariant obtained using the new Stokes solution and interpolated to markers and 
eii_stress is the stress-based second invariant of strain rate defined as follows:

    exx_stress = sxx2_marker/2.0/eta + dsxx_marker/2.0/dt/mu
    exy_stress = sxy2_marker/2.0/eta + dsxy_marker/2.0/dt/mu
    eii_stress = sqrt(exx_stress^2.0 + exy_stress^2.0)               eq. 2

where sxx2_marker and sxy2_marker are grid forecasted normal and shear stress 
interpolated to markers, dsxx_marker and dsxy_marker are the grid stress changes 
interpolated to markers, dt is the model time step, eta is viscosity and mu is 
the elastic modulus.

Equation 2 is derived from the visco-elastic stress forecast equations:

    sij = 2*eta*eij*(1-f) + sij_o*f                                    eq. 3

where sij is the stress, sij_o is the initial stress, eij is the strain rate, eta 
is the viscosity, f is a visco-elastic factor defined by the following equation:

    f = eta/(eta + mu*dt)                                              eq. 4

where mu is the shear modulus and dt is the model time step. The change in stress 
dsij can be obtained by subtracting the initial stress from equation 3:

    dsij = sij - sij_o = 2*eta*eii(1-f) - sij_o*(1-f),
    dsij = (2*eta*eii - sij_o)*(1-f).                                  eq. 5

Rearranging equation 5 gives the following expression for strain rate:

    eii = dsij/(2*eta) + sij_o/(2*eta) + dsij/(2*mu*dt)                eq. 6

Equation 6 can be reformulated using current stress sij as follows:

    eii = sij/(2*eta) + dsij/(2*mu*dt).                                eq. 7

See pgs 190 and 191 of Gerya (2010) for more details.

Updated Arrays
==============
model.markers.arrays.strain
--------------------------
marker_sr_ratio: Array{Float64,1}
    Ratio of nodal-stress-based to nodal-velocity-based strain rate invariant
    for each marker.
"""
function calculate_marker_strain_rate_ratio!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    timestep = model.timestep.parameters.main_time_loop.timestep.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array

    gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array
    ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array
    gridx_vy = model.grids.arrays.staggered_vy.gridx_vy.array
    xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array

    sxx2 = model.stokes_continuity.arrays.stress.sxx2.array
    sxy2 = model.stokes_continuity.arrays.stress.sxy2.array
    dsxx = model.stokes_continuity.arrays.stress_change.dsxx.array
    dsxy = model.stokes_continuity.arrays.stress_change.dsxy.array

    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_exx = model.markers.arrays.strain.marker_exx.array
    marker_exy = model.markers.arrays.strain.marker_exy.array
    marker_sr_ratio = model.markers.arrays.strain.marker_sr_ratio.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    mat_mu = model.materials.arrays.mat_mu.array

    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            # Define indexes and normalized distances for upper left pressure
            # grid node in the cell that contains the marker. Notes that
            # indices and distances for the basic grid have been pre-computed.
            @inbounds begin
                x_marker = marker_x[imarker] # marker x-location
                y_marker = marker_y[imarker] # marker y-location
                yn       = marker_yn[imarker] # upper left basic grid node y-index
                xn       = marker_xn[imarker] # upper left basic grid node x-index
                dy       = marker_dy[imarker] # normalized y-distance between marker and upper left basic grid node
                dx       = marker_dx[imarker] # normalized x-distance between marker and upper left basic grid node
            end
            yn_ul_pr = MarkerGridMapping.upr_left_index_pressure(ynum, yn, y_marker, gridy_vx)
            xn_ul_pr = MarkerGridMapping.upr_left_index_pressure(xnum, xn, x_marker, gridx_vy)
            dy_ul_pr = MarkerGridMapping.upr_left_dist_pressure(yn_ul_pr, y_marker, gridy_vx, ystp_vx)
            dx_ul_pr = MarkerGridMapping.upr_left_dist_pressure(xn_ul_pr, x_marker, gridx_vy, xstp_vy)
            # Interpolate normal grid stress to marker
            sxx2_marker = GridToMarker.get_marker_value(yn_ul_pr, xn_ul_pr, dy_ul_pr, dx_ul_pr, sxx2)
            # Interpolate shear grid stress to marker
            sxy_marker = GridToMarker.get_marker_value(yn, xn, dy, dx, sxy2)
            # Interpolate grid normal stress change to marker
            dsxx2_marker = GridToMarker.get_marker_value(yn_ul_pr, xn_ul_pr, dy_ul_pr, dx_ul_pr, dsxx)
            # Interpolate grid shear stress change to marker
            dsxy_marker = GridToMarker.get_marker_value(yn, xn, dy, dx, dsxy)
            # Get additional marker information
            @inbounds begin
                material_id    = marker_matid[imarker]
                shear_modulus  = mat_mu[material_id]
                viscosity      = marker_eta[imarker]
                strain_rate_xx = marker_exx[imarker]
                strain_rate_xy = marker_exy[imarker]
            end
            # Calculate the second strain rate invariant for the marker using
            # values interpolated from grid
            eii_velocity = second_invariant_strain_rate(strain_rate_xx, strain_rate_xy)
            # Calculate strain-rate ratio
            marker_sr_ratio[imarker] = marker_strain_rate_ratio(
                timestep, sxx2_marker, sxy_marker, dsxx2_marker, dsxy_marker,
                shear_modulus, viscosity, eii_velocity
            )
        end
    end
    return nothing
end

@inline function second_invariant_strain_rate(
    strain_rate_xx::Float64,
    strain_rate_xy::Float64
)::Float64
    eii = sqrt(strain_rate_xx*strain_rate_xx + strain_rate_xy*strain_rate_xy)
    return eii
end

""" Calculate marker strain rate ratio.

Returns
-------
marker_sr_ratio: Float64
    Ratio of nodal-stress-based to nodal-velocity-based strain rate invariant
    for each marker.
"""
function marker_strain_rate_ratio(
    timestep::Float64,
    sxx_marker::Float64,
    sxy_marker::Float64,
    dsxx_marker::Float64,
    dsxy_marker::Float64,
    shear_modulus::Float64,
    viscosity::Float64,
    eii_velocity::Float64
)::Float64
    if eii_velocity > 0
        if viscosity == 0
            @warn "viscosity is zero: $viscosity"
        end
        if timestep == 0
            @warn "timestep is zero: $timestep"
        end
        if shear_modulus == 0
            @warn "shear_modulus is zero: $shear_modulus"
        end
        if eii_velocity == 0
            @warn "eii_velocity is zero: $eii_velocity"
        end
        timestep_factor = 1.0/(timestep*shear_modulus)
        eii_stress = sqrt(
            (0.5*(sxx_marker/viscosity + dsxx_marker*timestep_factor))^2 +
            (0.5*(sxy_marker/viscosity + dsxy_marker*timestep_factor))^2
        )
        marker_sr_ratio = eii_stress/eii_velocity
    else
        marker_sr_ratio = 1.0
    end
    return marker_sr_ratio
end

end # module 