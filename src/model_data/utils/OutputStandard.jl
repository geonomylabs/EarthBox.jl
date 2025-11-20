module OutputStandard

using ..MarkerContainer: Markers
using ..HeatEquationContainer: HeatEquation
using ..StokesContinuityContainer: StokesContinuity

"""
    OutputLists

A struct that manages lists of marker and transport array objects that can be saved to file.
"""
mutable struct OutputLists
    marker_array_obj_list::Union{Nothing, Vector}
    transport_array_obj_list::Union{Nothing, Vector}

    function OutputLists(
        markers::Markers,
        heat_equation::HeatEquation,
        stokes_continuity::StokesContinuity
    )
        new_obj = new(nothing, nothing)
        make_marker_array_obj_output_list!(new_obj, markers)
        make_transport_array_obj_output_list!(new_obj, heat_equation, stokes_continuity)
        return new_obj
    end
end

"""
    make_marker_array_obj_output_list!(output_lists::OutputLists, markers::Markers)

Define list of marker array objects that may be saved to file.

This method generates a list of marker array objects that may be saved
to file if the output configuration is set to true.
"""
function make_marker_array_obj_output_list!(output_lists::OutputLists, markers::Markers)
    output_lists.marker_array_obj_list = [
        markers.arrays.location.marker_x,
        markers.arrays.location.marker_y,
        markers.arrays.thermal.marker_TK,
        markers.arrays.material.marker_matid,
        markers.arrays.material.marker_serpentinization,
        markers.arrays.material.marker_rho,
        markers.arrays.grid_marker_relationship.marker_xn,
        markers.arrays.grid_marker_relationship.marker_yn,
        markers.arrays.stress.marker_sxx,
        markers.arrays.stress.marker_sxy,
        markers.arrays.rheology.marker_eta,
        markers.arrays.strain.marker_exx,
        markers.arrays.strain.marker_exy,
        markers.arrays.pressure.marker_pr,
        markers.arrays.strain.marker_GII,
        markers.arrays.strain.marker_strain_plastic,
        markers.arrays.strain.marker_strain_rate_plastic,
        markers.arrays.strain.marker_sr_ratio,
        markers.arrays.melt.marker_meltfrac,
        markers.arrays.melt.marker_extracted_meltfrac,
        markers.arrays.melt.marker_extractable_meltfrac,
        markers.arrays.strat.marker_age,
        markers.arrays.rheology.marker_fric_ini,
        markers.arrays.rheology.marker_fric,
        markers.arrays.rheology.marker_cohesion,
        markers.arrays.rheology.marker_preexp,
        markers.arrays.rheology.marker_pfailure,
    ]
end

"""
    make_transport_array_obj_output_list!(output_lists::OutputLists, heat_equation::HeatEquation, stokes_continuity::StokesContinuityEquations)

Define list of transport array objects that will be saved to file.
"""
function make_transport_array_obj_output_list!(
    output_lists::OutputLists,
    heat_equation::HeatEquation,
    stokes_continuity::StokesContinuity
)
    output_lists.transport_array_obj_list = [
        heat_equation.arrays.temperature.tk1,
        stokes_continuity.arrays.viscosity.etas1,
        stokes_continuity.arrays.density.rho1,
        heat_equation.arrays.thermal_conductivity.kt1,
        stokes_continuity.arrays.stress.sxx2,
        stokes_continuity.arrays.stress.sxy2,
        stokes_continuity.arrays.pressure.pr1,
        stokes_continuity.arrays.strain_rate_and_spin.eii,
        stokes_continuity.arrays.plastic_def.plastics,
        stokes_continuity.arrays.plastic_def.plasticn
    ]
end

end # module OutputStandard 