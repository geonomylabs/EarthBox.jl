module Advection

include("options/Options.jl")
include("location/MarkersLocation.jl")
include("strain/MarkersStrain.jl")
include("stress_rotation/MarkersStressRotation.jl")
include("velocity_and_spin/VelocityAndSpin.jl")

using InteractiveUtils
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox: InitializationTools
import EarthBox: OptionTools
import EarthBox.Arrays: ArrayUtils
import .Options: OptionState
import .Options: get_options
import .Options: option_ids
import .Options: option_names

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_local_adaptive_time_stepping::Symbol
    marker_cell_displ_max::Symbol
    subgrid_diff_coef_temp::Symbol
    subgrid_diff_coef_stress::Symbol
end

function make_marker_advection_schemes_string()::String
    marker_advection_schemes_string = ""
    for (option_id, option_state) in get_options()
        option_name = Symbol(option_state.option_name)
        marker_advection_schemes_string *= """
## $(option_state.option_name)
- `advection_scheme` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
- **Description**: $(option_state.description)

"""
    end
    return marker_advection_schemes_string
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize marker advection scheme parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing the 
    model parameters and arrays.

# Keyword Arguments
- `advection_scheme::Union{Int, String, Nothing}=nothing`: Controls the type of marker advection scheme.
    See the **Marker Advection Schemes** section below for information on available marker advection schemes.
    The marker advection scheme is stored in the model data container as an integer ID (`itype_move_markers`) 
    and a corresponding string name (`stype_move_markers`). If `advection_scheme` is nothing the current 
    marker advection scheme defined in the model data container will be used. The marker advection scheme 
    parameters can be accessed from the model data container as follows:
    - `itype_move_markers = model.markers.parameters.advection.itype_move_markers.value`
    - `stype_move_markers = model.markers.parameters.advection.stype_move_markers.value`

# Keyword Arguments
- `$(PDATA.iuse_local_adaptive_time_stepping.name)::Int64`:
    - $(PDATA.iuse_local_adaptive_time_stepping.description)
- `$(PDATA.marker_cell_displ_max.name)::Float64`:
    - $(PDATA.marker_cell_displ_max.description)
- `$(PDATA.subgrid_diff_coef_temp.name)::Float64`:
    - $(PDATA.subgrid_diff_coef_temp.description)
- `$(PDATA.subgrid_diff_coef_stress.name)::Float64`:
    - $(PDATA.subgrid_diff_coef_stress.description)

---
# Marker Advection Schemes
---
$(make_marker_advection_schemes_string())
"""
function initialize!(
    model::ModelData;
    advection_scheme::Union{Int, String, Symbol, Nothing}=nothing,
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    _ = InitializationTools.update_option_id_using_input_option_name(
        get_options(), advection_scheme, model, get_option_id_from_model, update_option_id)
    return nothing
end

function get_max_runge_kutta_order(model::ModelData)::Int
    """ Define maximum Runge Kutta order using option id
    """
    return get_option_id_from_model(model)
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Marker Advection Option")
end

""" Advect markers and rotate stress tensor.

This function uses a 4-th order Runge Kutta scheme to interpolate velocity and spin to the 
marker locations. The marker locations and stress rotation are then updated using the 
interpolated values.
"""
function advect_markers_and_rotate_stress_tensor!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    markers_spin, markers_vx, markers_vy = interpolate_grid_velocity_and_spin_to_markers(model, inside_flags)
    update_marker_location!(model, inside_flags, markers_vx, markers_vy)
    update_marker_stress_rotation!(model, inside_flags, markers_spin)
    return nothing
end

""" Interpolate velocity and spin using a 4-th order Runge Kutta scheme.
"""
function interpolate_grid_velocity_and_spin_to_markers(model::ModelData, inside_flags::Vector{Int8})
    markers_spin = nothing
    markers_vx = nothing
    markers_vy = nothing
    max_runge_kutta_order = get_max_runge_kutta_order(model)
    if max_runge_kutta_order > 0
        @timeit_memit "Finished interpolating velocity and spin using Runge Kutta" begin
            (
                markers_vx, markers_vy, markers_spin
            ) = VelocityAndSpin.interpolate_using_runge_kutta(
                model, max_runge_kutta_order, inside_flags)
        end
    end
    return markers_spin, markers_vx, markers_vy
end

""" Update marker location using interpolated velocity.
"""
function update_marker_location!(
    model::ModelData,
    inside_flags::Vector{Int8},
    markers_vx::Union{Vector{Float64}, Nothing}=nothing,
    markers_vy::Union{Vector{Float64}, Nothing}=nothing
)::Nothing
    max_runge_kutta_order = get_max_runge_kutta_order(model)
    if max_runge_kutta_order > 0 && !isnothing(markers_vx) && !isnothing(markers_vy)
        @timeit_memit "Finished calculating new marker location" begin
            MarkersLocation.update!(model, inside_flags, markers_vx, markers_vy)
        end
    end

    return nothing
end

""" Update marker stress rotation using interpolated spin.
"""
function update_marker_stress_rotation!(
    model::ModelData,
    inside_flags::Vector{Int8},
    markers_spin::Union{Vector{Float64}, Nothing}=nothing
)::Nothing

    max_runge_kutta_order = get_max_runge_kutta_order(model)
    if max_runge_kutta_order > 0 && !isnothing(markers_spin)
        @timeit_memit "Finished updating marker stress rotation" begin
            MarkersStressRotation.update!(model, inside_flags, markers_spin)
        end
    end

    return nothing
end

""" Update marker strain
"""
function update_marker_strain!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "Finished updating marker strain" begin
        MarkersStrain.update_total_strain!(model, inside_flags)
        MarkersStrain.update_plastic_strain!(model, inside_flags)
    end
    return nothing
end

function get_option_id_from_model(model::ModelData)::Int
    return model.markers.parameters.advection.itype_move_markers.value
end

function get_stype_from_model(model::ModelData)::String
    return model.markers.parameters.advection.stype_move_markers.value
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.markers.parameters.advection.itype_move_markers.value = option_id
    return nothing
end

end # module 