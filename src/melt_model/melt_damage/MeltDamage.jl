module MeltDamage

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit, print_info
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.Markers.MarkerMaterials.MaterialGroupIDs: get_sticky_material_ids
import EarthBox.MathTools: zero_or_one
import ..MoltenZone: get_x_limits_of_molten_zone

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_melt_damage::Symbol
    iuse_probabilistic_melt_damage::Symbol
    melt_damage_distance::Symbol
    melt_damage_factor::Symbol
    melt_damage_taper_distance::Symbol
    maximum_damage_probability::Symbol
    magmatic_crust_height_threshold::Symbol
    magmatic_crust_height_minimum::Symbol
    magmatic_crust_height_maximum::Symbol
    magmatic_crust_height_intermediate::Symbol
    intermediate_damage_probability::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize melt damage model parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing the 
    model parameters and arrays.

# Keyword Arguments
- `$(PDATA.iuse_melt_damage.name)::Int64`
    - $(PDATA.iuse_melt_damage.description)
- `$(PDATA.iuse_probabilistic_melt_damage.name)::Int64`
    - $(PDATA.iuse_probabilistic_melt_damage.description)
- `$(PDATA.melt_damage_distance.name)::Float64`
    - $(PDATA.melt_damage_distance.description)
- `$(PDATA.melt_damage_factor.name)::Float64`
    - $(PDATA.melt_damage_factor.description)
- `$(PDATA.melt_damage_taper_distance.name)::Float64`
    - $(PDATA.melt_damage_taper_distance.description)
- `$(PDATA.maximum_damage_probability.name)::Float64`
    - $(PDATA.maximum_damage_probability.description)
- `$(PDATA.magmatic_crust_height_threshold.name)::Float64`
    - $(PDATA.magmatic_crust_height_threshold.description)
- `$(PDATA.magmatic_crust_height_minimum.name)::Float64`
    - $(PDATA.magmatic_crust_height_minimum.description)
- `$(PDATA.magmatic_crust_height_maximum.name)::Float64`
    - $(PDATA.magmatic_crust_height_maximum.description)
- `$(PDATA.magmatic_crust_height_intermediate.name)::Float64`
    - $(PDATA.magmatic_crust_height_intermediate.description)
- `$(PDATA.intermediate_damage_probability.name)::Float64`
    - $(PDATA.intermediate_damage_probability.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

function update_melt_damage!(model::ModelData)::Nothing
    iuse_melt_damage = model.materials.parameters.melt_damage.iuse_melt_damage.value
    if iuse_melt_damage == 1
        @timeit_memit "Finished updating markers for melt damage" begin
            reset_melt_damage!(model)
            if use_melt_damage_model(model) === true
                apply_melt_damage_to_drainage_basins!(model)
            end
            print_min_max_melt_damage(model)
        end
    end
    return nothing
end

function apply_melt_damage_to_drainage_basins!(model::ModelData)::Nothing
    print_info("Applying melt damage to drainage basins", level=1)
    avg_shallow_partial_melt_xcoors = 
        model.melting.arrays.extraction.avg_shallow_partial_melt_xcoors.array
    avg_shallow_partial_melt_ycoors = 
        model.melting.arrays.extraction.avg_shallow_partial_melt_ycoors.array
    iuse_probabilistic_melt_damage = 
        model.materials.parameters.melt_damage.iuse_probabilistic_melt_damage.value

    ndrainage_basin = model.melting.parameters.extraction.ndrainage_basin.value
    for idrainage_basin in 1:ndrainage_basin
        xmin_molten_domain, xmax_molten_domain = 
            get_x_limits_of_molten_zone(model, idrainage_basin)

        avg_shallow_partial_melt_xcoor_mantle = 
            avg_shallow_partial_melt_xcoors[idrainage_basin]
        avg_shallow_partial_melt_ycoor_mantle = 
            avg_shallow_partial_melt_ycoors[idrainage_basin]

        use_damage_marker_loop = use_melt_damage_loop(xmin_molten_domain, xmax_molten_domain)

        print_melt_damage_info(
            model, idrainage_basin, avg_shallow_partial_melt_xcoor_mantle,
            avg_shallow_partial_melt_ycoor_mantle, xmin_molten_domain,
            xmax_molten_domain, use_damage_marker_loop)

        if use_damage_marker_loop === true
            if iuse_probabilistic_melt_damage == 0
                melt_damage_marker_loop!(
                    model, avg_shallow_partial_melt_xcoor_mantle,
                    avg_shallow_partial_melt_ycoor_mantle)
            else
                melt_damage_marker_loop_probabilistic!(
                    model, avg_shallow_partial_melt_xcoor_mantle,
                    avg_shallow_partial_melt_ycoor_mantle)
            end
        end
    end
    return nothing
end

"""
Reset melt damage for all markers. A value of 1.0 indicates no damage since the 
friction coefficient is divided by this factor. Reset to 1 is only done if the 
melt damage model is activated.
"""
function reset_melt_damage!(model::ModelData)::Nothing
    iuse_melt_damage = 
        model.materials.parameters.melt_damage.iuse_melt_damage.value
    marker_melt_damage = model.markers.arrays.strain.marker_melt_damage.array
    marknum = model.markers.parameters.distribution.marknum.value
    if iuse_melt_damage == 1
        Threads.@threads for imarker in 1:marknum
            @inbounds marker_melt_damage[imarker] = 1.0
        end
    end
    return nothing
end

function melt_damage_marker_loop!(
    model::ModelData,
    avg_shallow_partial_melt_xcoor_mantle::Float64,
    avg_shallow_partial_melt_ycoor_mantle::Float64
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    melt_damage_distance = model.materials.parameters.melt_damage.melt_damage_distance.value
    melt_damage_factor = model.materials.parameters.melt_damage.melt_damage_factor.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    matids = model.markers.arrays.material.marker_matid.array
    marker_melt_damage = model.markers.arrays.strain.marker_melt_damage.array

    sticky_ids = get_sticky_material_ids(model)

    Threads.@threads for imarker in 1:marknum
        @inbounds begin
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            matid = matids[imarker]
        end
        is_sticky = check_sticky(matid, sticky_ids)
        in_damage_zone = check_damage_zone(
            x_marker, y_marker, melt_damage_distance,
            avg_shallow_partial_melt_xcoor_mantle,
            avg_shallow_partial_melt_ycoor_mantle
            )
        if is_sticky === false && in_damage_zone === true
            damage_factor = calculate_damage_factor_cos(
                x_marker, avg_shallow_partial_melt_xcoor_mantle,
                melt_damage_distance, melt_damage_factor
                )
            @inbounds marker_melt_damage[imarker] = damage_factor
        end
    end
    return nothing
end

""" Calculate probabilistic melt damage for each marker.

# Updated Arrays
- `marker_melt_damage`: Melt damage factor for each marker.
"""
function melt_damage_marker_loop_probabilistic!(
    model::ModelData,
    avg_shallow_partial_melt_xcoor_mantle::Float64,
    avg_shallow_partial_melt_ycoor_mantle::Float64
)::Nothing
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    matids = model.markers.arrays.material.marker_matid.array

    melt_damage_distance = 
        model.materials.parameters.melt_damage.melt_damage_distance.value
    melt_damage_factor = 
        model.materials.parameters.melt_damage.melt_damage_factor.value

    central_damage_probability = calculate_variable_melt_damage_probability(model)
    marker_melt_damage = model.markers.arrays.strain.marker_melt_damage.array
    sticky_ids = get_sticky_material_ids(model)

    marknum = model.markers.parameters.distribution.marknum.value
    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = matids[imarker]
        is_sticky = check_sticky(matid, sticky_ids)
        in_damage_zone = check_damage_zone(
            x_marker, y_marker, melt_damage_distance,
            avg_shallow_partial_melt_xcoor_mantle,
            avg_shallow_partial_melt_ycoor_mantle)
        if is_sticky === false && in_damage_zone === true
            damage_factor = calculate_damage_factor_probabilistic(
                x_marker, avg_shallow_partial_melt_xcoor_mantle,
                melt_damage_distance, melt_damage_factor,
                central_damage_probability)
            marker_melt_damage[imarker] = damage_factor
        end
    end
    return nothing
end

function calculate_variable_melt_damage_probability(model::ModelData)::Float64
    characteristic_magmatic_crust_height = 
        model.melting.parameters.extrusion.characteristic_magmatic_crust_height.value
    magmatic_crust_height_threshold = 
        model.materials.parameters.melt_damage.magmatic_crust_height_threshold.value
    magmatic_crust_height_minimum = 
        model.materials.parameters.melt_damage.magmatic_crust_height_minimum.value
    magmatic_crust_height_maximum = 
        model.materials.parameters.melt_damage.magmatic_crust_height_maximum.value
    magmatic_crust_height_intermediate = 
        model.materials.parameters.melt_damage.magmatic_crust_height_intermediate.value
    maximum_melt_damage_probability = 
        model.materials.parameters.melt_damage.maximum_damage_probability.value
    intermediate_damage_probability = 
        model.materials.parameters.melt_damage.intermediate_damage_probability.value

    melt_damage_probability = linear_melt_damage_probability_model(
        characteristic_magmatic_crust_height, magmatic_crust_height_threshold,
        magmatic_crust_height_minimum, magmatic_crust_height_intermediate,
        magmatic_crust_height_maximum, maximum_melt_damage_probability,
        intermediate_damage_probability
        )

    print_variable_melt_damage_probability_info(
        characteristic_magmatic_crust_height, magmatic_crust_height_threshold,
        magmatic_crust_height_minimum, magmatic_crust_height_maximum,
        maximum_melt_damage_probability, melt_damage_probability)

    return melt_damage_probability
end

function calculate_variable_melt_damage_factor(model::ModelData)::Float64
    characteristic_magmatic_crust_height = 
        model.melting.parameters.extrusion.characteristic_magmatic_crust_height.value
    melt_damage_factor_maximum = 
        model.materials.parameters.melt_damage.melt_damage_factor.value
    magmatic_crust_height_threshold = 
        model.materials.parameters.melt_damage.magmatic_crust_height_threshold.value
    magmatic_crust_height_minimum = 
        model.materials.parameters.melt_damage.magmatic_crust_height_minimum.value
    magmatic_crust_height_maximum = 
        model.materials.parameters.melt_damage.magmatic_crust_height_maximum.value

    melt_damage_factor = linear_melt_damage_model(
        characteristic_magmatic_crust_height, magmatic_crust_height_threshold,
        magmatic_crust_height_minimum, magmatic_crust_height_maximum,
        melt_damage_factor_maximum
        )

    return melt_damage_factor
end

function print_variable_melt_damage_probability_info(
    characteristic_magmatic_crust_height::Float64,
    magmatic_crust_height_threshold::Float64,
    magmatic_crust_height_minimum::Float64,
    magmatic_crust_height_maximum::Float64,
    maximum_damage_probability::Float64,
    central_damage_probability::Float64
)::Nothing
    println(">> Variable melt damage probability info")
    println(">>    Characteristic magmatic crust height: ", 
            characteristic_magmatic_crust_height)
    println(">>    Maximum damage probability: ", maximum_damage_probability)
    println(">>    Magmatic crust height threshold: ", magmatic_crust_height_threshold)
    println(">>    Magmatic crust height minimum: ", magmatic_crust_height_minimum)
    println(">>    Magmatic crust height maximum: ", magmatic_crust_height_maximum)
    println(">>    Central damage probability: ", central_damage_probability)
    return nothing
end

function linear_melt_damage_model(
    characteristic_magmatic_crust_height::Float64,
    magmatic_crust_height_threshold::Float64,
    magmatic_crust_height_minimum::Float64,
    magmatic_crust_height_maximum::Float64,
    melt_damage_factor_maximum::Float64
)::Float64
    if characteristic_magmatic_crust_height < magmatic_crust_height_threshold
        melt_damage_factor = 1.0
    elseif characteristic_magmatic_crust_height <= magmatic_crust_height_minimum
        melt_damage_factor = 1.0
    elseif characteristic_magmatic_crust_height >= magmatic_crust_height_maximum
        melt_damage_factor = melt_damage_factor_maximum
    else
        melt_damage_factor = 1.0 + 
            (melt_damage_factor_maximum - 1.0) /
            (magmatic_crust_height_maximum - magmatic_crust_height_minimum) *
            (characteristic_magmatic_crust_height - magmatic_crust_height_minimum)
    end
    return melt_damage_factor
end

function linear_melt_damage_probability_model(
    characteristic_magmatic_crust_height::Float64,
    magmatic_crust_height_threshold::Float64,
    magmatic_crust_height_minimum::Float64,
    magmatic_crust_height_intermediate::Float64,
    magmatic_crust_height_maximum::Float64,
    maximum_damage_probability::Float64,
    intermediate_damage_probability::Float64
)::Float64
    if characteristic_magmatic_crust_height < magmatic_crust_height_threshold
        damage_probability = 0.0
    elseif characteristic_magmatic_crust_height <= magmatic_crust_height_minimum
        damage_probability = 0.0
    elseif characteristic_magmatic_crust_height >= magmatic_crust_height_maximum
        damage_probability = maximum_damage_probability
    elseif characteristic_magmatic_crust_height <= magmatic_crust_height_intermediate
        damage_probability = intermediate_damage_probability /
            (magmatic_crust_height_intermediate - magmatic_crust_height_minimum) *
            (characteristic_magmatic_crust_height - magmatic_crust_height_minimum)
    else
        damage_probability = intermediate_damage_probability +
            (maximum_damage_probability - intermediate_damage_probability) /
            (magmatic_crust_height_maximum - magmatic_crust_height_intermediate) *
            (characteristic_magmatic_crust_height - magmatic_crust_height_intermediate)
    end
    return damage_probability
end

@inline function check_sticky(matid::Int16, sticky_ids::Tuple{Int16, Int16})::Bool
    return matid == sticky_ids[1] || matid == sticky_ids[2]
end

function print_melt_damage_info(
    model::ModelData,
    idrainage_basin::Int,
    avg_shallow_partial_melt_xcoor_mantle::Float64,
    avg_shallow_partial_melt_ycoor_mantle::Float64,
    xmin_molten_domain::Float64,
    xmax_molten_domain::Float64,
    use_damage_marker_loop::Bool
)::Nothing
    characteristic_magmatic_crust_height = 
        model.melting.parameters.extrusion.characteristic_magmatic_crust_height.value
    iuse_melt_damage = 
        model.materials.parameters.melt_damage.iuse_melt_damage.value

    print_info("Melt damage info for drainage basin: $idrainage_basin", level=2)
    print_info("iuse_melt_damage: $iuse_melt_damage", level=3)
    print_info("Characteristic magmatic crust height: $characteristic_magmatic_crust_height", level=3)
    print_info("Avg shallow partial melt x-coor mantle: $avg_shallow_partial_melt_xcoor_mantle", level=3)
    print_info("Avg shallow partial melt y-coor mantle: $avg_shallow_partial_melt_ycoor_mantle", level=3)
    print_info("X-min molten domain: $xmin_molten_domain", level=3)
    print_info("X-max molten domain: $xmax_molten_domain", level=3)
    print_info("Use melt damage marker loop for basin: $use_damage_marker_loop", level=3)
    return nothing
end

function print_min_max_melt_damage(model::ModelData)::Nothing
    marker_melt_damage = model.markers.arrays.strain.marker_melt_damage.array
    min_damage = minimum(marker_melt_damage)
    max_damage = maximum(marker_melt_damage)
    print_info("Min melt damage: $min_damage", level=2)
    print_info("Max melt damage: $max_damage", level=2)
    return nothing
end

@inline function calculate_damage_factor_cos(
    x_marker::Float64,
    avg_shallow_partial_melt_xcoor_mantle::Float64,
    melt_damage_distance::Float64,
    melt_damage_factor::Float64
)::Float64
    xmin = avg_shallow_partial_melt_xcoor_mantle - melt_damage_distance
    xmax = avg_shallow_partial_melt_xcoor_mantle + melt_damage_distance
    damage_width = xmax - xmin
    if x_marker < xmin || x_marker > xmax
        damage_factor = 1.0
    else
        theta_prime = (x_marker - avg_shallow_partial_melt_xcoor_mantle) / damage_width * 2π
        damage_factor = 0.5*(melt_damage_factor - 1.0) * (cos(theta_prime) + 1.0) + 1.0
    end
    return damage_factor
end

@inline function check_damage_zone(
    x_marker::Float64,
    y_marker::Float64,
    melt_damage_distance::Float64,
    x_shallow_partially_molten::Float64,
    y_shallow_partially_molten::Float64
)::Bool
    xmin = x_shallow_partially_molten - melt_damage_distance
    xmax = x_shallow_partially_molten + melt_damage_distance
    return xmin < x_marker < xmax && y_marker < y_shallow_partially_molten
end

""" Check if melt damage model should be used to modify properties.

The model is only used if the x-coordinates of the gabbroic molten zone are
positive (i.e. gabbroic molten material is present in the model).

# Arguments
- `xmin_molten_domain`: Minimum x-coordinate of molten domain
- `xmax_molten_domain`: Maximum x-coordinate of molten domain

# Returns
- `check`: Flag indicating whether melt damage loop should be used to modify 
  properties
"""
function use_melt_damage_loop(
    xmin_molten_domain::Float64,
    xmax_molten_domain::Float64
)::Bool
    return !(xmin_molten_domain < 0.0 && xmax_molten_domain < 0.0)
end

""" Check if melt damage model should be used to modify properties.

The model is only used if the melt damage option is activated and the 
characteristic magmatic crust height is greater than a threshold.

# Arguments
- `model`: The model data container

# Returns
- `check`: Flag indicating whether melt damage loop should be used to modify 
  properties
"""
function use_melt_damage_model(model::ModelData)::Bool
    characteristic_magmatic_crust_height = 
        model.melting.parameters.extrusion.characteristic_magmatic_crust_height.value
    magmatic_crust_height_threshold = 
        model.materials.parameters.melt_damage.magmatic_crust_height_threshold.value

    iuse_melt_damage = model.materials.parameters.melt_damage.iuse_melt_damage.value
    check = false
    if iuse_melt_damage == 1
        check = true
    end
    if characteristic_magmatic_crust_height < magmatic_crust_height_threshold
        check = false
    end
    if check === true
        print_info("Using melt damage model", level=2)
    else
        print_info("Not using melt damage model", level=2)
    end
    return check
end

function calculate_damage_factor_probabilistic(
    x_marker::Float64,
    avg_shallow_partial_melt_xcoor_mantle::Float64,
    melt_damage_distance::Float64,
    melt_damage_factor::Float64,
    central_damage_probability::Float64
)::Float64
    damage_probability = calculate_damage_probability(
        x_marker, avg_shallow_partial_melt_xcoor_mantle,
        melt_damage_distance, central_damage_probability
        )
    result = zero_or_one(damage_probability)
    return result == 1.0 ? melt_damage_factor : 1.0
end

""" Calculate spatial factor for marker melt damage.

# Arguments
- `x_marker`: x-coordinate of marker
- `avg_shallow_partial_melt_xcoor_mantle`: x-coordinate of shallow partially 
  molten mantle
- `melt_damage_distance`: Melt damage distance
- `central_damage_probability`: Maximum probability of damage in the central 
  part of the damage zone

# Returns
- `damage_probability`: Damage probability
"""
function calculate_damage_probability(
    x_marker::Float64,
    avg_shallow_partial_melt_xcoor_mantle::Float64,
    melt_damage_distance::Float64,
    central_damage_probability::Float64
)::Float64
    xmin = avg_shallow_partial_melt_xcoor_mantle - melt_damage_distance
    xmax = avg_shallow_partial_melt_xcoor_mantle + melt_damage_distance
    damage_width = xmax - xmin
    if x_marker < xmin || x_marker > xmax
        damage_probability = 0.0
    else
        theta_prime = (x_marker - avg_shallow_partial_melt_xcoor_mantle) /
            damage_width * 2π
        damage_probability = central_damage_probability / 2.0 * (cos(theta_prime) + 1.0)
    end
    return damage_probability
end

end # module 