module XdmfFields

include("core/TimeStep.jl")
include("core/TimeSteps.jl")

using JLD2
import JLD2: attributes
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer.OutputStandard: OutputLists
import ..OutputDTypes: Mesh2djld
import ..OutputDTypes: ScalarField
import ..XdmfUtils: get_xdmf_number_type_for_array
import ..XdmfUtils: getoutform
import .TimeSteps: FieldsXdmfTimeSteps
import .TimeStep: FieldsXdmfTimeStep

function export_xdmf_fields(
    fields_xdmf_steps::FieldsXdmfTimeSteps,
    model::ModelData,
    output_lists::OutputLists,
    gridx_b_km::Vector{Float64},
    gridy_b_km::Vector{Float64},
    gridx_pr_km::Vector{Float64},
    gridy_pr_km::Vector{Float64}
)::Nothing
    noutput = model.timestep.parameters.output_steps.noutput.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_Myr = model.conversion.parameters.sec_per_Myr.value

    scalars_on_nodes = ScalarField[]
    scalars_on_cells = ScalarField[]
    ifind_basic_grid = false
    mesh2djld = nothing

    for obj in output_lists.transport_array_obj_list
        if obj.grid_type in ["basic"] && !ifind_basic_grid
            gridy_km = copy(gridy_b_km)
            gridx_km = copy(gridx_b_km)

            ynum = length(gridy_km)
            xnum = length(gridx_km)
            jld_mesh_filename = "fields_$(lpad(noutput, 5, "0")).jld"

            mesh_data = Dict(
                "ynum" => ynum,
                "xnum" => xnum,
                "ysize" => gridy_km[end],
                "xsize" => gridx_km[end],
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
            ifind_basic_grid = true
        end

        if obj.grid_type in ["basic"]
            push!(scalars_on_nodes, ScalarField(
                obj.outform.fname,
                obj.outform.fname,
                get_xdmf_number_type_for_array(obj.array),
                obj.outform.units,
                "None",
                "basic",
                getoutform(obj)
            ))
        elseif obj.grid_type in ["pressure"]
            push!(scalars_on_cells, ScalarField(
                obj.outform.fname,
                obj.outform.fname,
                get_xdmf_number_type_for_array(obj.array),
                obj.outform.units,
                "None",
                "pressure",
                getoutform(obj)
            ))
        end
    end

    model_time = timesum/sec_per_Myr
    time_units = "Myr"
    fields_xdmf_step = FieldsXdmfTimeStep(
        model_time, time_units, noutput,
        mesh2djld, scalars_on_nodes, scalars_on_cells
    )
    TimeStep.make_jld2_file(fields_xdmf_step, fields_xdmf_steps.output_dir)
    string = TimeStep.get_xdmf_string_for_timestep(fields_xdmf_step)
    TimeSteps.add_step_to_xdmf_string!(fields_xdmf_steps, string)
    TimeSteps.save(fields_xdmf_steps)
    return nothing
end

end # module 