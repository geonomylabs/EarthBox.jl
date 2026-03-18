module CentralWeakening

using Random
import EarthBox.ModelDataContainer: ModelData
import EarthBox.MathTools: zero_or_one
import EarthBox.StaggeredGrid: check_for_ttype_refinement
import ...FrictionRandomizer: weaken_initial_friction_coefficient
import EarthBox.DebugPlots: plot_marker_fric

"""
    initialize!(model)

Update marker friction coefficients via central weakening.

# Updated Arrays
## Updated arrays from group `model.markers.arrays.rheology`
- `marker_friction.array::Vector{Float64}`:
    - Array of friction coefficients (sine of friction angle) of markers.
- `marker_fric_ini.array::Vector{Float64}`:
    - Array of initial friction coefficients (sine of friction angle) of markers.
"""
function initialize!(model::ModelData)
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    marker_random = rand(mxnum*mynum)

    marker_x = model.markers.arrays.location.marker_x.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_plastic = model.materials.arrays.mat_plastic.array

    random_friction = model.materials.parameters.random_friction
    central_delta_fric_coef = random_friction.central_delta_fric_coef.value
    central_weakening_probability = random_friction.central_weakening_probability.value
    
    marker_fric_ini = model.markers.arrays.rheology.marker_fric_ini.array
    marker_fric = model.markers.arrays.rheology.marker_fric.array

    matids_sticky_air = model.materials.dicts.matid_types["StickyAir"]
    matids_sticky_water = model.materials.dicts.matid_types["StickyWater"]
    
    xsize = model.grids.parameters.geometry.xsize.value
    if check_for_ttype_refinement(model)
        xo_highres = model.grids.parameters.refinement.xo_highres.value
        xf_highres = model.grids.parameters.refinement.xf_highres.value
        if !isnan(xo_highres) && !isnan(xf_highres)
            xmin = xo_highres
            xmax = xf_highres
        else
            xmin = 0.0
            xmax = xsize
        end
    else
        xmin = 0.0
        xmax = xsize
    end

    println("CentralWeakening: xmin = $xmin, xmax = $xmax")

    marknum = model.markers.parameters.distribution.marknum.value
    for imarker in 1:marknum
        matid = marker_matid[imarker]
        x_marker = marker_x[imarker]
        friction_coefficient = mat_plastic[matid, 3]
        # Avoid sticky air and water
        if matid in matids_sticky_air || matid in matids_sticky_water
            marker_fric_ini[imarker] = friction_coefficient
        else
            if x_marker >= xmin && x_marker <= xmax
                weakening_probability = calculate_weakening_probability(
                    x_marker, xmin, xmax, central_weakening_probability)
                result = zero_or_one(weakening_probability)
                if result == 1.0
                    random_number = marker_random[imarker]
                    friction_coefficient = weaken_initial_friction_coefficient(
                        friction_coefficient, central_delta_fric_coef, random_number)
                end
            end
            marker_fric_ini[imarker] = friction_coefficient

        end
        marker_fric[imarker] = marker_fric_ini[imarker]
    end

    iplot_marker_fric = 0
    if iplot_marker_fric == 1
        plot_marker_fric(model)
    end

    return nothing
end


""" Weakening probability as a function of distance from zone center.

Probability is maximum at the center of the zone [xmin, xmax] and smoothly
decreases to zero at the edges via a half-cosine profile.

# Arguments
- `x_marker`: x-coordinate of marker
- `xmin`: Left bound of the weakening zone
- `xmax`: Right bound of the weakening zone
- `central_weakening_probability`: Maximum weakening probability at the center

# Returns
- `weakening_probability`: Weakening probability (0 outside [xmin, xmax])
"""
function calculate_weakening_probability(
    x_marker::Float64,
    xmin::Float64,
    xmax::Float64,
    central_weakening_probability::Float64
)::Float64
    damage_width = xmax - xmin
    center = (xmin + xmax) / 2.0
    if x_marker < xmin || x_marker > xmax
        weakening_probability = 0.0
    else
        theta_prime = (x_marker - center) / damage_width * 2π
        weakening_probability = central_weakening_probability / 2.0 * (cos(theta_prime) + 1.0)
    end
    return weakening_probability
end

end # module 