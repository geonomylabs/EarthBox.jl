module XdmfTopography

include("core/GetTopoData.jl")
include("core/TimeStep.jl")
include("core/TimeSteps.jl")

import EarthBox.ModelDataContainer: ModelData
import ..OutputDTypes: Topo2djld
import .GetTopoData: get_topo_data
import .TimeSteps: TopographyXdmfTimeSteps

function export_topography(steps::TopographyXdmfTimeSteps, model::ModelData)
    iuse_topo = model.topography.parameters.model_options.iuse_topo.value
    if iuse_topo == 1
        topo_data = get_topo_data(model)
        topo2djld = Topo2djld(topo_data)
        
        model_time = topo_data["time"]
        time_units = topo_data["time_units"]
        noutput = topo_data["noutput"]
        
        topography_xdmf_time_step = TimeStep.TopographyXdmfTimeStep(
            topo2djld, model_time, time_units, noutput)
            
        TimeStep.make_jld2_file(topography_xdmf_time_step, steps.output_dir)
        string = TimeStep.get_xdmf_string_for_timestep(topography_xdmf_time_step)
        TimeSteps.add_step_to_xdmf_string(steps, string)
        TimeSteps.save(steps)
    end
end

end # module 