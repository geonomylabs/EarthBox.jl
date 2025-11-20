module Stopping

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

struct ValidParameterNames
    iuse_boundary_displacement::Symbol
    displ_limit::Symbol
    iuse_extensional_strain::Symbol
    strain_limit::Symbol
end

function initialize(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidParameterNames); kwargs...)
    return nothing
end

""" Test for cases where a break from the time loop is required.

Breaking is associated with cases with boundary displacement limits
and cases with extensional strain limits.

Returns
-------
break_check::Bool
    Boolean denoting cases where a break in the time loop is required.
"""
function terminate_loop(model::ModelData)::Bool
    break_check = false
    boundary_dis = model.timestep.parameters.boundary_displacement_stopping
    if boundary_dis.iuse_boundary_displacement.value == 1
        break_check = check_boundary_displacement_stopping(model, break_check)
    end
    if boundary_dis.iuse_extensional_strain.value == 1
        break_check = check_extensional_strain_stopping(model, break_check)
    end
    return break_check
end

""" Check boundary displacement stopping criteria.
"""
function check_boundary_displacement_stopping(
    model::ModelData,
    break_check::Bool
)::Bool
    boundary_dis = model.timestep.parameters.boundary_displacement_stopping
    displ_limit = boundary_dis.displ_limit.value
    dxshort = calc_boundary_displacement(model)
    boundary_displacement_msg(dxshort, displ_limit)
    if dxshort > displ_limit
        boundary_displacement_limit_msg()
        break_check = true
    end
    return break_check
end

""" Check extensional strain stopping criteria.
"""
function check_extensional_strain_stopping(
    model::ModelData,
    break_check::Bool
)::Bool
    boundary_dis = model.timestep.parameters.boundary_displacement_stopping
    strain_limit = boundary_dis.strain_limit.value
    displacement_total, strain_total = calc_extensional_displacement_and_strain(model)
    extensional_strain_msg(model, strain_total, displacement_total)
    
    if strain_total >= strain_limit
        extensional_strain_limit_msg(model, strain_total)
        extensional_strain_limit_reached_msg()
        break_check = true
    end
    
    return break_check
end

""" Print extensional strain message.
"""
function extensional_strain_msg(
    model::ModelData,
    strain_total::Float64,
    displacement_total::Float64
)
    boundary_dis = model.timestep.parameters.boundary_displacement_stopping
    strain_limit = boundary_dis.strain_limit.value
    
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    
    println(
        ">> Extensional strain: $strain_total : Limit: $strain_limit ",
        "Displacement_total (m): $displacement_total ",
        "ntimestep: $ntimestep : Time (Myr): $(timesum/sec_per_myr)"
    )
end

"""
    extensional_strain_limit_reached_msg()

Print extensional strain limit reached message.
"""
function extensional_strain_limit_reached_msg()
    println(">> Extensional strain limit reached. Breaking from time loop.")
end

"""
    extensional_strain_limit_msg(
        pymodel::ModelData,
        strain_total::Float64
    )

Print strain limit.
"""
function extensional_strain_limit_msg(
    model::ModelData,
    strain_total::Float64
)
    boundary_dis = model.timestep.parameters.boundary_displacement_stopping
    strain_limit = boundary_dis.strain_limit.value
    
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    
    println(
        ">> Total strain $strain_total reached ",
        "strain limit $strain_limit. ",
        "Breaking from time loop at $(timesum/sec_per_myr) Myr."
    )
end

function get_strain_limit(model::ModelData)::Float64
    boundary_dis = model.timestep.parameters.boundary_displacement_stopping
    return boundary_dis.strain_limit.value
end

function boundary_displacement_msg(dxshort::Float64, displ_limit::Float64)
    println(
        ">> Total boundary displacement (dxshort) (m): $dxshort",
        " limit (m): $displ_limit"
    )
end

function boundary_displacement_limit_msg()
    println(">> Boundary displacement limit reached. Breaking from time loop.")
end

""" Calculate and print boundary displacement.

Returns
-------
dxshort::Float64
    Boundary displacement in meters.
"""
function calc_boundary_displacement(model::ModelData)::Float64
    velocity_internal_x = model.bcs.parameters.velocity.velocity_internal_x.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    dxshort = abs(velocity_internal_x * timesum)
    return dxshort
end

""" Calculate and print extensional displacement and strain.

Returns
-------
displacement_total::Float64
    Total extensional displacement in meters.
strain_total::Float64
    Total extensional strain.
"""
function calc_extensional_displacement_and_strain(
    model::ModelData
)::Tuple{Float64, Float64}
    xsize_start = model.grids.parameters.geometry.xsize_start.value
    full_velocity_extension = model.bcs.parameters.velocity.full_velocity_extension.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    
    displacement_total = full_velocity_extension * timesum
    strain_total = displacement_total / xsize_start
    
    return displacement_total, strain_total
end

end # module 