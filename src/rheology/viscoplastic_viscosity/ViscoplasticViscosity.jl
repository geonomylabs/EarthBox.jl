module ViscoplasticViscosity

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: check_in_domain

"""
    update_marker_viscoplastic_viscosity!(model::ModelData)::Nothing

This function updates marker viscoplastic viscosity.

If using a node-based plasticity, marker viscosity is only updated at 
initialization or when marker viscosity is zero. If using a marker-based 
plasticity, marker viscosity is always updated.

# Arguments
- `model`: Main model container object.

# Updated Arrays
- `model.markers.arrays.rheology.marker_eta.array`: Marker visco-plastic viscosity (Pa.s)
"""
function update_marker_viscoplastic_viscosity!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    itype_global = model.stokes_continuity.parameters.picard.itype_global.value
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    marknum = model.markers.parameters.distribution.marknum.value
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array
    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                viscosity_flow = marker_eta_flow[imarker]
                viscosity_viscoplastic = marker_eta[imarker]
            end
            @inbounds marker_eta[imarker] = get_updated_marker_viscosity(
                ntimestep, itype_global, viscosity_flow, viscosity_viscoplastic)
        end
    end
    return nothing
end

@inline function get_updated_marker_viscosity(
    ntimestep::Int,
    itype_global::Int,
    viscosity_flow::Float64,
    viscosity_viscoplastic::Float64
)::Float64
    # Always update marker viscosity if using marker-based plastic global loop
    if itype_global == 0
        viscosity_viscoplastic = viscosity_flow
    # Only update marker viscosity with flow viscosity at initialization when
    # using node-base plastic global loop viscosity
    elseif itype_global == 1
        if ntimestep == 0 || viscosity_viscoplastic == 0.0
            viscosity_viscoplastic = viscosity_flow
        end
    end
    return viscosity_viscoplastic
end

end # module 