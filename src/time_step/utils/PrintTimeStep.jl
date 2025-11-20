module PrintTimeStep

using EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: print_info

function print_timestep(model::ModelData, msg::String)
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    print_info("Current time step (Myr) $(msg) : $(round(timestep/sec_per_myr, digits=6))", level=2)
end

end # module 