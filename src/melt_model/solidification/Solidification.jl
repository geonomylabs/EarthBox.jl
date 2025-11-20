module Solidification

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Markers.MarkerFriction.FrictionRandomizer: randomize_initial_friction_coefficient
import EarthBox.SurfaceProcesses: calculate_age_ma

""" Transform magma and lava markers to solid if below the liquidus.

This function loops over markers and calls functions that update the marker
material id array ``marker_matid`` to account for solidification of
molten material when melt fraction falls below 1. These functions also
reset several marker arrays to account for the change in composition.

# Returns
- `Float64`: Maximum depth of solidified crust above region of melting (meters)
"""
function solidify!(model::ModelData)::Float64
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array

    matid_types = model.materials.dicts.matid_types
    # Add checks or warnings if more than 1 type is defined for the following:
    if length(matid_types["SolidifiedGabbro"]) > 0
        matid_gabbro = matid_types["SolidifiedGabbro"][1]
    else
        matid_gabbro = -1
    end
    if length(matid_types["ExtractedGabbroicMagma"]) > 0
        matid_gabbroic_magma = matid_types["ExtractedGabbroicMagma"][1]
    else
        matid_gabbroic_magma = -1
    end
    if length(matid_types["SolidifiedLayeredGabbro"]) > 0
        matid_layered_gabbro = matid_types["SolidifiedLayeredGabbro"][1]
    else
        matid_layered_gabbro = -1
    end
    if length(matid_types["ExtractedLayeredGabbroicMagma"]) > 0
        matid_layered_gabbroic_magma = matid_types["ExtractedLayeredGabbroicMagma"][1]
    else
        matid_layered_gabbroic_magma = -1
    end
    if length(matid_types["SolidifiedBasalt"]) > 0
        matid_basalt = matid_types["SolidifiedBasalt"][1]
    else
        matid_basalt = -1
    end
    if length(matid_types["ExtrudedGabbroicMagma"]) > 0
        matid_lava = matid_types["ExtrudedGabbroicMagma"][1]
    else
        matid_lava = -1
    end
    age_ma = calculate_age_ma(model)

    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    marker_random = rand(mxnum * mynum)

    iuse_random_fric = model.materials.parameters.random_friction.iuse_random_fric.value
    delta_fric_coef = model.materials.parameters.random_friction.delta_fric_coef.value

    mat_plastic = model.materials.arrays.mat_plastic.array

    marknum = model.markers.parameters.distribution.marknum.value
    Threads.@threads for imarker in 1:marknum
        @inbounds begin
            matid = marker_matid[imarker]
            meltfrac = marker_meltfrac[imarker]
        end
        # Gabbroic magma to solidified gabbro
        transform_molten_material_to_solid(
            model, imarker, matid, matid_gabbroic_magma, matid_gabbro, meltfrac,
            marker_random, mat_plastic, iuse_random_fric, delta_fric_coef,
            age_ma
        )
        # Layered gabbroic magma to solidified layered gabbro
        transform_molten_material_to_solid(
            model, imarker, matid, matid_layered_gabbroic_magma,
            matid_layered_gabbro, meltfrac, marker_random, mat_plastic,
            iuse_random_fric, delta_fric_coef, age_ma
        )
        # Lava to solidified basalt
        transform_molten_material_to_solid(
            model, imarker, matid, matid_lava, matid_basalt, meltfrac,
            marker_random, mat_plastic, iuse_random_fric, delta_fric_coef,
            age_ma
        )
    end

    # This is currently not being used
    ymax_solidified_crust = max_depth_of_solidified_crust_above_melting_region(model)

    return ymax_solidified_crust
end

""" Get updated friction coefficient.

# Arguments
- `mat_plastic::Matrix{Float64}`: Material plastic properties array
- `matid::Int16`: Material ID
- `iuse_random_fric::Int64`: Whether to use random friction
- `delta_fric_coef::Float64`: Delta friction coefficient
- `random_number::Float64`: Random number for friction calculation

# Returns
- `Float64`: Updated friction coefficient
"""
@inline function get_updated_friction_coefficient(
    mat_plastic::Matrix{Float64},
    matid::Int16,
    iuse_random_fric::Int64,
    delta_fric_coef::Float64,
    random_number::Float64
)::Float64
    # Get friction coefficient from material model
    friction_coefficient = mat_plastic[matid, 3]
    if iuse_random_fric == 1
        friction_coefficient = randomize_initial_friction_coefficient(
            friction_coefficient, delta_fric_coef, random_number)
    end
    return friction_coefficient
end

""" Transform molten material to solidified material.

# Arguments
- `model::ModelData`: Model data container object
- `imarker::Int64`: Marker index
- `matid::Int16`: Material ID
- `matid_molten::Int16`: Molten material ID
- `matid_solid::Int16`: Solid material ID
- `meltfrac::Float64`: Melt fraction
- `marker_random::Vector{Float64}`: Random numbers for markers
- `mat_plastic::Matrix{Float64}`: Material plastic properties array
- `iuse_random_fric::Int64`: Whether to use random friction
- `delta_fric_coef::Float64`: Delta friction coefficient
- `age_ma::Float64`: Age in million years
"""
@inline function transform_molten_material_to_solid(
    model::ModelData,
    imarker::Int64,
    matid::Int16,
    matid_molten::Int16,
    matid_solid::Int16,
    meltfrac::Float64,
    marker_random::Vector{Float64},
    mat_plastic::Matrix{Float64},
    iuse_random_fric::Int64,
    delta_fric_coef::Float64,
    age_ma::Float64
)::Nothing
    if below_liquidus(matid, matid_molten, meltfrac)
        random_number = marker_random[imarker]
        friction_coefficient = get_updated_friction_coefficient(
            mat_plastic, matid_solid, iuse_random_fric, delta_fric_coef,
            random_number
        )
        transform_marker_to_solid(
            model, imarker, matid_solid, age_ma, friction_coefficient)
    end
    return nothing
end

""" Check for extruded molten mantle material below liquidus.

# Arguments
- `matid::Int16`: Material ID
- `matid_molten::Int16`: Molten material ID
- `meltfrac::Float64`: Melt fraction

# Returns
- `Bool`: True if extruded mantle is below liquidus
"""
@inline function below_liquidus(
    matid::Int16,
    matid_molten::Int16,
    meltfrac::Float64
)::Bool
    return matid == matid_molten && meltfrac < 1
end

""" Transform magma marker to solidified material.

# Updated Arrays
- `model.markers.arrays.stress.marker_sxx`: Normal stress of markers in Pascals
- `model.markers.arrays.stress.marker_sxy`: Shear stress of markers in Pascals
- `model.markers.arrays.rheology.marker_eta`: Effective viscosity of markers in Pa.s
- `model.markers.arrays.strain.marker_exx`: Normal strain rate of markers in 1/s
- `model.markers.arrays.strain.marker_exy`: Shear strain rate of markers in 1/s
- `model.markers.arrays.strain.marker_GII`: Strain of markers
- `model.markers.arrays.strain.marker_strain_plastic`: Plastic strain of markers
- `model.markers.arrays.strain.marker_sr_ratio`: Ratio of marker strain rate 
  calculated using stress change and Maxwell model over strain rate interpolated 
  from grid
- `model.materials.arrays.marker_matid`: Material ID of markers
- `model.markers.arrays.melt.marker_meltfrac`: Melt fraction of markers
- `model.markers.arrays.strat.marker_age`: Age of formation of markers in Ma

# Arguments
- `model::ModelData`: Model data container object
- `imarker::Int64`: Marker index
- `matid_solid::Int16`: Solid material ID
- `age_ma::Float64`: Age in million years
- `friction_coeff::Float64`: Friction coefficient
"""
@inline function transform_marker_to_solid(
    model::ModelData,
    imarker::Int64,
    matid_solid::Int16,
    age_ma::Float64,
    friction_coeff::Float64
)::Nothing
    @inbounds begin
        # set to material ID to a material ID for solidified material
        model.markers.arrays.material.marker_matid.array[imarker] = Int16(matid_solid)
        # Reset strain rate ratio
        model.markers.arrays.strain.marker_sr_ratio.array[imarker] = 1
        # Reset strain
        model.markers.arrays.strain.marker_GII.array[imarker] = 0
        model.markers.arrays.strain.marker_strain_plastic.array[imarker] = 0
        # Reset stress
        model.markers.arrays.stress.marker_sxx.array[imarker] = 0.0
        model.markers.arrays.stress.marker_sxy.array[imarker] = 0.0
        # Reset viscosity
        model.markers.arrays.rheology.marker_eta.array[imarker] = 0.0
        # Reset strain rate
        model.markers.arrays.strain.marker_exx.array[imarker] = 0.0
        model.markers.arrays.strain.marker_exy.array[imarker] = 0.0
        # Reset friction coefficient
        model.markers.arrays.rheology.marker_fric_ini.array[imarker] = friction_coeff
        model.markers.arrays.rheology.marker_fric.array[imarker] = friction_coeff
        # Define age of melting
        model.markers.arrays.strat.marker_age.array[imarker] = age_ma
    end
    return nothing
end

""" Determine the max depth of solidified crust above region of melting.

# Arguments
- `model::ModelData`: Model data container object

# Returns
- `Float64`: Maximum depth of solidified crust above region of melting (meters)
"""
function max_depth_of_solidified_crust_above_melting_region(
    model::ModelData
)::Float64
    marknum = model.markers.parameters.distribution.marknum.value
    xmid_mol = model.melting.parameters.extraction.xmid_mol.value
    width_mol = model.melting.parameters.extraction.width_mol.value
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    matid_types = model.materials.dicts.matid_types

    if length(matid_types["SolidifiedGabbro"]) > 0
        matid_gabbro = matid_types["SolidifiedGabbro"][1]
    else
        matid_gabbro = -1
    end
    if length(matid_types["SolidifiedBasalt"]) > 0
        matid_basalt = matid_types["SolidifiedBasalt"][1]
    else
        matid_basalt = -1
    end

    ymax_solid = -1e32
    xwidth_check = width_mol * 2
    xmin_check = xmid_mol - xwidth_check / 2.0
    xmax_check = xmid_mol + xwidth_check / 2.0

    for imarker in 1:marknum
        @inbounds begin
            matid = marker_matid[imarker]
            y_marker = marker_y[imarker]
            x_marker = marker_x[imarker]
        end
        if solidified_crust(
            matid, x_marker, y_marker, xmin_check, xmax_check, 
            ymax_solid, matid_gabbro, matid_basalt
        )
            ymax_solid = y_marker
        end
    end
    return ymax_solid
end

@inline function solidified_crust(
    matid::Int16,
    x_marker::Float64,
    y_marker::Float64,
    xmin_check::Float64,
    xmax_check::Float64,
    ymax_solid::Float64,
    matid_gabbro::Int16,
    matid_basalt::Int16,
)::Bool
    return matid in (matid_gabbro, matid_basalt) && xmin_check < x_marker < xmax_check && y_marker > ymax_solid
end

end # module 