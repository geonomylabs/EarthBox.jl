module ModelParametersInfo

include("FormatValueStrings.jl")

import DataStructures: OrderedDict
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParameterRegistry: get_input_parameter_name_list
import .FormatValueStrings: get_value_string_for_terminal

function print_parameters(model)
    object_list = [
        model,
        model.grids.parameters.geometry,
        model.grids.parameters.grid_options,
        model.grids.parameters.refinement,
        model.markers.parameters.advection,
        model.markers.parameters.distribution,
        model.markers.parameters.recycling,
        model.markers.parameters.subgrid_diffusion,
        model.timestep.parameters.boundary_displacement_stopping,
        model.bcs.parameters.bc_options,
        model.bcs.parameters.pressure,
        model.bcs.parameters.temperature,
        model.bcs.parameters.velocity,
        model.bcs.parameters.velocity_step,
        model.stokes_continuity.parameters.build,
        model.stokes_continuity.parameters,
        model.stokes_continuity.parameters.output,
        model.stokes_continuity.parameters.picard,
        model.stokes_continuity.parameters.residual_norms,
        model.stokes_continuity.parameters.solution_norms,
        model.stokes_continuity.parameters.velocity_calc_option,
        model.materials.parameters.material_description,
        model.materials.parameters.random_friction,
        model.materials.parameters.softening,
        model.materials.parameters.stress_limits_power_law,
        model.materials.parameters.stress_limits_yield,
        model.materials.parameters.viscosity_limits,
        model.materials.parameters.compaction,
        model.materials.parameters.hydrothermal,
        model.geometry.parameters.crustal_hole,
        model.geometry.parameters.earth_layering,
        model.geometry.parameters.litho_strong_zones,
        model.geometry.parameters.plume,
        model.geometry.parameters.rayleigh_taylor,
        model.geometry.parameters.sticky_air_geometry,
        model.geometry.parameters.weak_fault,
        model.geometry.parameters.weak_seed,
        model.heat_equation.parameters,
        model.melting.parameters.extraction,
        model.melting.parameters.extrusion,
        model.melting.parameters.melting_material_ids,
        model.melting.parameters.options,
        model.melting.parameters.rheology,
        model.topography.parameters.depo_and_erosion_rates,
        model.topography.parameters.downhill_diffusion,
        model.topography.parameters.model_options,
        model.topography.parameters.sealevel,
        model.topography.parameters.topo_grid,
        model.interpolation.parameters.interp_options,
        model.gravity.parameters,
        model.carbonate.parameters,
        model.benchmarks.parameters
    ]

    master_input_names, _ = get_input_parameter_name_list()

    param_dict = make_parameter_info_dict(
        object_list, 1, master_input_names)

   print_param_dict(param_dict)
end

function make_parameter_info_dict(
    members_list::Vector{Any},
    iuse_input_name_filter::Int,
    master_input_names::Vector{String}
)::OrderedDict{String, Vector{Any}}
    param_dict = OrderedDict{String, Vector{Any}}()
    for member in members_list
        for attr_name in propertynames(member)
            attr = getproperty(member, attr_name)
            collect_parameter_info(
                attr_name,
                attr,
                iuse_input_name_filter,
                master_input_names,
                param_dict
            )
        end
    end
    return param_dict
end

function collect_parameter_info(
    attr_name::Symbol,
    attr::Any,
    iuse_input_name_filter::Int,
    master_input_names::Vector{String},
    param_dict::OrderedDict{String, Vector{Any}}
)::Nothing
    if attr isa Union{ParameterFloat, ParameterInt, ParameterStr}
        name = String(attr_name)
        iuse_name = 1
        if iuse_input_name_filter == 1
            if name in master_input_names
                iuse_name = 1
            else
                iuse_name = 0
            end
        end
        if iuse_name == 1
            value_str = get_value_string_for_terminal(attr.value)
            param_dict[name] = [
                value_str,
                attr.name,
                attr.units,
                attr.description
            ]
        end
    end
    return nothing
end

function print_param_dict(param_dict::OrderedDict{String, Vector{Any}})
    print_info("", level=1)
    print_info("Model Parameters", level=1)
    print_info("----------------", level=1)
    for (key, value) in param_dict
        print_info("$(rpad(key, 45)) : $(lpad(value[1], 12)) : $(rpad(value[3], 8)) : $(value[4])", level=2)
    end
    print_info("", level=1)
end

end # module 