"""
    XdmfTimeStepsManager

Module responsible for managing XDMF output files and timesteps.

This module handles:
- Writing model data to XDMF format files
- Managing separate XDMF files for fields, markers, velocity and topography
- Converting and formatting data for output
- Coordinating timestep information across output files

The XdmfTimeStepsManager acts as the central coordinator for XDMF output,
organizing the writing of different data types to their respective files
while maintaining consistent timestep information.
"""
module XdmfTimeStepsManager

include("utils/OutputDTypes.jl")
include("utils/XdmfParts.jl")
include("utils/XdmfUtils.jl")
include("utils/XdmfTimeStepTools.jl")
include("fields/XdmfFields.jl")
include("markers/XdmfMarkers.jl")
include("velocity/XdmfVelocity.jl")
include("topography/XdmfTopography.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConfigurationManager.OutputConfig: OutputConfigState
import EarthBox.ModelDataContainer.OutputStandard: OutputLists
import .XdmfFields.TimeSteps: FieldsXdmfTimeSteps
import .XdmfMarkers.TimeSteps: MarkersXdmfTimeSteps
import .XdmfVelocity.TimeSteps: VelocityXdmfTimeSteps
import .XdmfTopography.TimeSteps: TopographyXdmfTimeSteps

struct XdmfTimeSteps
    fields_xdmf_steps::FieldsXdmfTimeSteps
    markers_xdmf_steps::MarkersXdmfTimeSteps
    velocity_xdmf_steps::VelocityXdmfTimeSteps
    topography_xdmf_steps::TopographyXdmfTimeSteps
end

function XdmfTimeSteps(output_dir::String)
    return XdmfTimeSteps(
        FieldsXdmfTimeSteps(output_dir),
        MarkersXdmfTimeSteps(output_dir),
        VelocityXdmfTimeSteps(output_dir),
        TopographyXdmfTimeSteps(output_dir)
    )
end

function convert_grid_arrays_to_km_for_output(
    model::ModelData
)::Tuple{
    Vector{Float64},
    Vector{Float64},
    Vector{Float64},
    Vector{Float64},
    Vector{Float64},
    Vector{Float64}
}
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    garrays = model.grids.arrays
    gridx_b_km = copy(garrays.basic.gridx_b.array) ./ 1000.0
    gridy_b_km = copy(garrays.basic.gridy_b.array) ./ 1000.0
    gridx_vy_km = copy(garrays.staggered_vy.gridx_vy.array[2:xnum]) ./ 1000.0
    gridy_vx_km = copy(garrays.staggered_vx.gridy_vx.array[2:ynum]) ./ 1000.0
    gridx_pr_km = copy(garrays.pressure.gridx_pr.array) ./ 1000.0
    gridy_pr_km = copy(garrays.pressure.gridy_pr.array) ./ 1000.0
    
    return gridx_b_km, gridy_b_km, gridx_vy_km, gridy_vx_km, gridx_pr_km, gridy_pr_km
end

function export_xdmf(
    xdmf_steps::XdmfTimeSteps,
    model::ModelData,
    output_lists::OutputLists,
    output_config::OutputConfigState
)
    XdmfMarkers.export_xdmf_markers(
        xdmf_steps.markers_xdmf_steps, model, output_lists, output_config)
    
    (
        gridx_b_km, gridy_b_km, gridx_vy_km, gridy_vx_km, gridx_pr_km, gridy_pr_km
    ) = convert_grid_arrays_to_km_for_output(model)
    
    XdmfFields.export_xdmf_fields(
        xdmf_steps.fields_xdmf_steps, model, output_lists, 
        gridx_b_km, gridy_b_km, gridx_pr_km, gridy_pr_km)
    
    XdmfVelocity.export_xdmf_velocity(
        xdmf_steps.velocity_xdmf_steps, model, 
        gridx_b_km, gridy_b_km, gridx_pr_km, gridy_pr_km)
    
    XdmfTopography.export_topography(
        xdmf_steps.topography_xdmf_steps, model)
end

end # module 