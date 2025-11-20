module CheckMagmaFlush

import EarthBox.ModelDataContainer: ModelData

function check_for_initial_magma_flush(model::ModelData)::Bool
    check = false
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    initial_magma_flush_steps = 
        model.melting.parameters.extrusion.initial_magma_flush_steps.value
    
    if initial_magma_flush_steps > 0
        if ntimestep <= initial_magma_flush_steps
            check = true
        end
    end
    return check
end

end # module 