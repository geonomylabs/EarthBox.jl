module Rheology

include("composite_creep/CompositeCreep.jl")
include("viscous_strain_softening/ViscousStrainSoftening.jl")
include("flow_viscosity/FlowViscosity.jl")
include("viscoplastic_viscosity/ViscoplasticViscosity.jl")
include("plastic_failure/PlasticFailure.jl")

import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: ModelData
import .PlasticFailure: MarkerPlasticity
import .PlasticFailure.MarkerPlasticity: option_names as marker_plasticity_names
import .PlasticFailure: FailureProperties
import .FlowViscosity: FlowViscosity

""" Initialize Rheology

This currently only initializes the marker_plasticity model.
"""
function initialize!(
    model::ModelData;
    marker_plasticity_model::Union{Int,String,Symbol,Nothing} = nothing
)::Nothing
    MarkerPlasticity.initialize!(
        model, marker_plasticity_model=marker_plasticity_model)
    return nothing
end

""" Update marker rheology properties.

Update plastic failure properties (cohesion and friction coefficients) 
and flow viscosity for markers including the effects of serpentinization 
if applicable. Marker arrays `marker_cohesion`, `marker_fric`,
`marker_eta_flow`, and `marker_preexp` are updated.

"""
function update_marker_rheology!(
    model::ModelData,
    boolean_options::Dict{String,Bool},
    inside_flags::Vector{Int8}
    #serpentinization_model::SerpentinizationModel
)::Nothing
    @timeit_memit "Finished calculating cohesion and friction coefficients" begin
        (
            use_boundary_friction_plasticity_model_sandbox
        ) = get(boolean_options, "use_boundary_friction_plasticity_model_sandbox", false)
        FailureProperties.calculate_failure_properties!(
            model, inside_flags, use_boundary_friction_plasticity_model_sandbox)
    end
    @timeit_memit "Finished updating marker flow viscosity" begin
        FlowViscosity.update_marker_flow_viscosity!(model, inside_flags)
    end
    # To Do: Add serpentinization model
    #update_marker_failure_properties(serpentinization_model, model)
    #update_marker_flow_viscosity(serpentinization_model, model)
end

""" Update marker viscoplastic viscosity

Update marker viscoplastic viscosity array `marker_eta`. This
array is initialized using the flow viscosity array
`marker_eta_flow` and updated for partial melt.

If marker-based plasticity us used, the viscoplastic viscosity is 
updated for yielding. If the viscoplastic stress forecast option is
selected (instead of elastic stress forecast) with marker-based 
plasticity the marker stress arrays `marker_sxx` and `marker_syy` 
are updated for yielding. 
"""
function update_marker_viscoplastic_viscosity!(
    model::ModelData,
    boolean_options::Dict{String,Bool},
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "Finished initial update of marker viscoplastic viscosity" begin
        ViscoplasticViscosity.update_marker_viscoplastic_viscosity!(model, inside_flags)
    end
    update_marker_viscoplastic_viscosity_for_yielding(model, boolean_options, inside_flags)
end

""" Update marker viscoplastic viscosity for yielding.

This method is only applied for marker-based plasticity as opposed
to grid-based plasticity.

Update marker viscoplastic viscosity array (marker_eta) and marker
plastic yielding flag (marker_pfailure) for plastic failure if
marker-based global loop is used (itype_plasticity = 0). Marker stress
(marker_sxx, marker_sxy) is also updated for marker-based global loop
if the viscoelastic marker plasticity option is selected as opposed
to the pure elastic option.
"""
function update_marker_viscoplastic_viscosity_for_yielding(
    model::ModelData,
    boolean_options::Dict{String,Bool},
    inside_flags::Vector{Int8}
)::Nothing
    no_yielding_in_mobile_wall = get(boolean_options, "no_yielding_in_mobile_wall", false)
    use_marker_plasticity = get(boolean_options, "use_marker_plasticity", false)
    if use_marker_plasticity
        @timeit_memit "Finished updating marker viscoplastic viscosity for yielding" begin
            MarkerPlasticity.update_plastic_yielding!(
                model, inside_flags, no_yielding_in_mobile_wall)
        end
    end
    return nothing
end

end # module

