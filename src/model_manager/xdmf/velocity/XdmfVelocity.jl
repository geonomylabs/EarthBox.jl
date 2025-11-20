module XdmfVelocity

include("core/TimeStep.jl")
include("core/TimeSteps.jl")

using JLD2
import EarthBox.ModelDataContainer: ModelData
import EarthBox.EarthBoxDtypes: AbstractXDMFTimeStep
import ..OutputDTypes: Mesh2djld, Vector2D
import ..XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ..XdmfParts: get_xdmf_end_for_collection_grid_and_domain
import ..XdmfParts: get_xdmf_start_of_timestep_grid
import ..XdmfParts: get_xdmf_end_of_timestep_grid
import ..XdmfParts: get_xdmf_topology_2drectmesh
import ..XdmfParts: get_xdmf_geometry_2d_vxvy
import ..XdmfParts: get_xdmf_vector_attribute_on_nodes
import ..XdmfUtils: getoutform
import ..XdmfUtils: get_xdmf_number_type_for_array
import ..XdmfTimeStepTools
import .TimeSteps: VelocityXdmfTimeSteps

function export_xdmf_velocity(
    steps::VelocityXdmfTimeSteps,
    model::ModelData,
    gridx_b_km::Array{Float64},
    gridy_b_km::Array{Float64},
    gridx_pr_km::Array{Float64},
    gridy_pr_km::Array{Float64}
)::Nothing
    noutput = model.timestep.parameters.output_steps.noutput.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    basic_vel = model.stokes_continuity.arrays.basic_grid_velocity
    units = basic_vel.vxb.outform.units

    gridy_km = copy(gridy_b_km)
    gridx_km = copy(gridx_b_km)

    ynum = size(gridy_km, 1)
    xnum = size(gridx_km, 1)
    jld_mesh_filename = "vel_cmyr_$(lpad(noutput, 5, "0")).jld"

    mesh_data = Dict(
        "ynum" => ynum,
        "xnum" => xnum,
        "ysize" => gridy_km[ynum],
        "xsize" => gridx_km[xnum],
        "jld_mesh_filename" => jld_mesh_filename,
        "jld_dataname_1dy" => "gridy_km",
        "jld_dataname_1dx" => "gridx_km",
        "jld_dataname_pr1dy" => "gridy_pr_km",
        "jld_dataname_pr1dx" => "gridx_pr_km",
        "gridy_km_array" => gridy_km,
        "gridx_km_array" => gridx_km,
        "gridy_pr_km_array" => gridy_pr_km,
        "gridx_pr_km_array" => gridx_pr_km
    )

    mesh2djld = Mesh2djld(mesh_data)

    vectorx = getoutform(basic_vel.vxb)
    vectory = getoutform(basic_vel.vyb)

    velocity = Vector2D(
        "vel_cm_yr",
        "velxyz_cmyr",
        "velx_cmyr",
        "vely_cmyr",
        "velmag_cmyr",
        get_xdmf_number_type_for_array(vectorx),
        units,
        "None",
        "basic",
        vectorx,
        vectory
    )

    model_time = timesum / sec_per_myr
    time_units = "Myr"

    velocity_xdmf_time_step = TimeStep.VelocityXdmfTimeStep(
        model_time,
        time_units,
        noutput,
        mesh2djld,
        velocity
    )

    TimeStep.make_jld2_file(velocity_xdmf_time_step, steps.output_dir)
    string = TimeStep.get_xdmf_string_for_timestep(velocity_xdmf_time_step)
    TimeSteps.add_step_to_xdmf_string(steps, string)
    TimeSteps.save_xdmf_string(steps)
    
    return nothing
end

end # module 