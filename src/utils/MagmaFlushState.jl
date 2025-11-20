module MagmaFlushState

import EarthBox.ModelDataContainer: ModelData

function update_use_extrusion_for_magma_flush(
    model::ModelData,
    iuse_extrusion::Int64
)::Int64
    use_magma_flush = check_magma_flush_state(model)
    if use_magma_flush
        iuse_extrusion = 1
    end
    return iuse_extrusion
end

function check_magma_flush_state(
    model::ModelData
)::Bool
    initial_magma_flush_steps = model.melting.parameters.extrusion.initial_magma_flush_steps.value
    check = false
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    if initial_magma_flush_steps > 0
        if ntimestep <= initial_magma_flush_steps
            check = true
        end
    end
    return check
end

end # module 