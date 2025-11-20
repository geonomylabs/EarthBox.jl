module XdmfMarkers

include("core/GetMarkerData.jl")
include("core/TimeStep.jl")
include("core/TimeSteps.jl")

using JLD2
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer.OutputStandard: OutputLists
import EarthBox.ConfigurationManager.OutputConfig: OutputConfigState
import ..OutputDTypes: Markers2djld
import ..OutputDTypes: ScalarField
import .GetMarkerData: get_marker_data
import .TimeSteps: MarkersXdmfTimeSteps

function export_xdmf_markers(
    markers_xdmf::MarkersXdmfTimeSteps,
    model::ModelData,
    output_lists::OutputLists,
    output_config::OutputConfigState
)
    marker_data = get_marker_data(model, output_lists, output_config)
    markers2djld = Markers2djld(marker_data)
    scalars_on_markers = marker_data["scalars_on_markers"]

    model_time = marker_data["time"]
    time_units = marker_data["time_units"]
    noutput = marker_data["noutput"]
    y_sealevel = marker_data["y_sealevel"]
    base_level_shift = marker_data["base_level_shift"]

    markers_xdmf_time_step = TimeStep.MarkersXdmfTimeStep(
        markers2djld, model_time, time_units, noutput,
        y_sealevel, base_level_shift, scalars_on_markers
    )

    TimeStep.make_jld2_file(markers_xdmf_time_step, markers_xdmf.output_dir)
    string = TimeStep.get_xdmf_string_for_timestep(markers_xdmf_time_step)
    TimeSteps.add_step_to_xdmf_string(markers_xdmf, string)
    TimeSteps.save(markers_xdmf)
end

end # module 