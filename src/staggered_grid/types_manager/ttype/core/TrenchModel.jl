module TrenchModel

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import ..Initialization

function update_highres_region_using_trench!(model::ModelData)::Float64
    update_trench_location_for_adaptive_mesh_refinement!(model)
    trench_location_m = model.grids.parameters.refinement.trench_location.value
    xo_highres = model.grids.parameters.refinement.xo_highres.value
    xf_highres = model.grids.parameters.refinement.xf_highres.value
    width_highres = xf_highres - xo_highres
    return trench_location_m - width_highres / 2.0
end

function update_trench_location_for_adaptive_mesh_refinement!(model::ModelData)::Float64
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value

    y_marker_top = 15000.0
    y_marker_base = 25000.0  # Bottom of domain used to track trench location, original value = 20 km

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    x_trench_limit = calc_x_limit_for_trench(model)

    domains = model.materials.dicts.matid_domains
    matid_asthenosphere = domains["Asthenosphere"]
    matid_weak_mantle_lithosphere = domains["WeakMantleLithosphere"]

    fracture_zone = model.geometry.parameters.fracture_zone
    x_fracture_zone_end = fracture_zone.x_fracture_zone_end.value

    # Initialize trench location
    if ntimestep == 0
        trench_location = x_fracture_zone_end
    else
        trench_location = model.grids.parameters.refinement.trench_location.value
    end

    x_highres_min, x_highres_max = get_x_limits_for_highres_region(model)

    use_only_weak_lith = true

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        if in_trench_region(x_marker, y_marker, x_highres_min, x_highres_max, y_marker_top, y_marker_base)
            matid = marker_matid[imarker]
            if !use_only_weak_lith
                check_fz_mantle = is_fracture_zone_mantle(
                    y_marker, y_marker_base, matid, matid_asthenosphere, matid_weak_mantle_lithosphere)
            else
                check_fz_mantle = is_fracture_zone_mantle_weak_lithosphere_only(
                    y_marker, y_marker_top, y_marker_base, matid, matid_weak_mantle_lithosphere
                )
            end
            if check_fz_mantle && x_marker > trench_location
                # && x_marker <= x_trench_limit
                trench_location = x_marker
            end
        end
    end
    model.grids.parameters.refinement.trench_location.value = trench_location
    print_info("Updated trench location (m): $(trench_location)", level=2)

    return trench_location
end

function get_x_limits_for_highres_region(model::ModelData)::Tuple{Float64, Float64}
    ixo_highres, ixf_highres = get_indices_of_high_res_region(model)
    gridx_b = model.grids.arrays.basic.gridx_b.array
    x_highres_min = gridx_b[ixo_highres]
    x_highres_max = gridx_b[ixf_highres]
    return x_highres_min, x_highres_max
end

function in_trench_region(
    x_marker::Float64, 
    y_marker::Float64,
    x_highres_min::Float64,
    x_highres_max::Float64,
    y_marker_top::Float64, 
    y_marker_base::Float64
)::Bool
    if y_marker > y_marker_top && y_marker < y_marker_base && x_marker > x_highres_min && x_marker < x_highres_max
        return true
    else
        return false
    end
end

"""
    calc_x_limit_for_trench(model::ModelData)::Float64

Calculate the current x-limit for the trench.
This limit is just to the right of the mid-point of the high-resolution region.
"""
function calc_x_limit_for_trench(model::ModelData)::Float64
    ixo_highres, ixf_highres = get_indices_of_high_res_region(model)
    ix_trench_limit = ixo_highres + Int64(floor(0.5*(ixf_highres - ixo_highres)))
    gridx_b = model.grids.arrays.basic.gridx_b.array
    x_trench_limit = gridx_b[ix_trench_limit]
    return x_trench_limit
end

"""
    get_indices_of_high_res_region(model::ModelData)::Tuple{Int64,Int64}

Get indices of high resolution region.
"""
function get_indices_of_high_res_region(model::ModelData)::Tuple{Int64,Int64}
    xo_highres = model.grids.parameters.refinement.xo_highres.value
    xf_highres = model.grids.parameters.refinement.xf_highres.value
    dx_highres = model.grids.parameters.refinement.dx_highres.value
    xnum = model.grids.parameters.geometry.xnum.value

    ixo_highres = Initialization.get_left_edge_x_index_of_highres_zone(
        xnum, xo_highres, xf_highres, dx_highres)
    ixf_highres = Initialization.calculate_final_xgrid_index_of_highres(
        ixo_highres, xo_highres, xf_highres, dx_highres)

    return (ixo_highres, ixf_highres)
end

"""
    is_fracture_zone_mantle(y_marker::Float64, y_marker_base::Float64, matid::Int16,
                          matid_asthenosphere::Int16, matid_weak_mantle_lithosphere::Int16)::Bool

Check if material is fracture zone mantle close to surface.
"""
function is_fracture_zone_mantle(
    y_marker::Float64,
    y_marker_base::Float64,
    matid::Int16,
    matid_asthenosphere::Int16,
    matid_weak_mantle_lithosphere::Int16
)::Bool
    check = false
    if matid == matid_asthenosphere || matid == matid_weak_mantle_lithosphere
        check = true
    end
    if y_marker > y_marker_base
        check = false
    end
    return check
end

function is_fracture_zone_mantle_weak_lithosphere_only(
    y_marker::Float64,
    y_marker_top::Float64,
    y_marker_base::Float64,
    matid::Int16,
    matid_weak_mantle_lithosphere::Int16
)::Bool
    check = false
    if matid == matid_weak_mantle_lithosphere
        check = true
    end
    if !(y_marker_top < y_marker < y_marker_base)
        check = false
    end
    return check
end

end # module
