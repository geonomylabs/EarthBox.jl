module SaltDeposition

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: load_parameters!
import ..Topography.AdvectionTools: advect_topography_vertically_for_salt_deposition!
import ..MarkerTransformation
import ..TopoInterpolation

export initialize!

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_salt_deposition::Symbol
    salt_start_time::Symbol
    salt_end_time::Symbol
    salt_deposition_rate::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Apply salt deposition model.

# Arguments
- `model::`[`ModelData`](@ref ModelData)
    - The model data container containing the model parameters and arrays.

# Keyword Arguments
- `$(PDATA.iuse_salt_deposition.name)::Int64`
    - $(PDATA.iuse_salt_deposition.description)
- `$(PDATA.salt_start_time.name)::Float64`
    - $(PDATA.salt_start_time.description)
- `$(PDATA.salt_end_time.name)::Float64`
    - $(PDATA.salt_end_time.description)
- `$(PDATA.salt_deposition_rate.name)::Float64`
    - $(PDATA.salt_deposition_rate.description)

"""
function initialize!(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

function apply_salt_deposition!(model::ModelData)::Nothing
    iuse_topo = model.topography.parameters.model_options.iuse_topo.value
    iuse_salt_deposition = 
        model.topography.parameters.model_options.iuse_salt_deposition.value

    if iuse_topo == 1 && iuse_salt_deposition == 1
        if check_salt_parameters(model)
            @timeit_memit "Finished applying salt deposition model" begin
                advect_topography_vertically_for_salt_deposition!(model)
                apply_marker_transformation_for_salt_deposition!(model)
            end
        end
    end
    return nothing
end

function check_salt_parameters(model::ModelData)::Bool
    matids_salt = model.materials.dicts.matid_types["Salt"]
    salt_check = true
    
    if isempty(matids_salt)
        salt_check = false
        println("!!! WARNING!!! Salt material not found. ",
                "Salt deposition model will not be applied.")
    else
        if matids_salt[1] == -1
            salt_check = false
            println("!!! WARNING!!! Salt material not found. ",
                    "Salt deposition model will not be applied.")
        end
    end
    
    return salt_check
end 

function apply_marker_transformation_for_salt_deposition!(
    model::ModelData
)::Nothing
    age_ma = calculate_age_ma(model)

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_matid = model.markers.arrays.material.marker_matid.array

    y_sealevel = model.topography.parameters.sealevel.y_sealevel
    dx_topo = model.topography.parameters.topo_grid.dx_topo
    toponum = model.topography.parameters.topo_grid.toponum
    gridt = model.topography.arrays.gridt

    matids_sticky_air = model.materials.dicts.matid_types["StickyAir"]
    matids_sticky_water = model.materials.dicts.matid_types["StickyWater"]
    matids_salt = model.materials.dicts.matid_types["Salt"]

    marknum = model.markers.parameters.distribution.marknum
    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        # Only transform markers below sea level
        if y_marker > y_sealevel
            ixn, dx_dist = TopoInterpolation.find_closest_topography_node_to_marker(
                x_marker, gridt, dx_topo, toponum)

            y_topo = TopoInterpolation.calculate_topography(dx_dist, gridt, ixn)

            if air_in_sediment(
                marker_matid[imarker], y_marker, y_topo, matids_sticky_air
            )
                MarkerTransformation.transform_marker_to_salt!(
                    model, imarker, age_ma, matids_salt
                )
            end

            if water_in_sediment(
                marker_matid[imarker], y_marker, y_topo, matids_sticky_water
            )
                MarkerTransformation.transform_marker_to_salt!(
                    model, imarker, age_ma, matids_salt
                )
            end
        end
    end
end

end # module