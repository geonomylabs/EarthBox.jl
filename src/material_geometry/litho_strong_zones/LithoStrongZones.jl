module LithoStrongZones

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.StaggeredGrid: check_for_ttype_refinement

export initialize!

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_strong_zones::Symbol
    x_left_strong::Symbol
    x_right_strong::Symbol
end

""" Initialize lithosphere strong zones geometry.

This material geometry involves a central lithospheric zone with different frictional plastic
properties than the surrounding lithosphere. The left and right boundaries of the central 
lithospheric zone are controlled by the parameters `x_left_strong` and `x_right_strong`, respectively.

Material properties of the central and lateral zones are defined using the material domains for the 
upper continental crust strong zone, lower continental crust strong zone, and lithospheric mantle 
strong zone in either a `material.yaml` file or `Material.jl` script.

# Keyword arguments:
- `iuse_strong_zones::Int`:
    - $(PDATA.iuse_strong_zones.description)
- `x_left_strong::Union{Float64, Nothing}`:
    - $(PDATA.x_left_strong.description)
- `x_right_strong::Union{Float64, Nothing}`:
    - $(PDATA.x_right_strong.description)

"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    update_parameters!(model)
    return nothing
end

""" Update lithosphere strong zones parameters using t-type refinement parameters.

If t-type refinement is selected and `xo_highres` and `xf_highres` are set, update the
lithosphere strong zones parameters using the t-type refinement parameters.

"""
function update_parameters!(model::ModelData)
    if conditions_for_update(model)
        xo_highres = model.grids.parameters.refinement.xo_highres.value
        xf_highres = model.grids.parameters.refinement.xf_highres.value
        
        litho_strong_zones = model.geometry.parameters.litho_strong_zones
        litho_strong_zones.x_left_strong.value = xo_highres
        litho_strong_zones.x_right_strong.value = xf_highres

        print_info("Updated lithosphere strong zones parameters using t-type refinement parameters.", level=1)
        print_info("Inputs:", level=2)
        print_info("xo_highres: $(xo_highres)", level=3)
        print_info("xf_highres: $(xf_highres)", level=3)
        print_info("Outputs:", level=2)
        print_info("x_left_strong: $(litho_strong_zones.x_left_strong.value)", level=3)
        print_info("x_right_strong: $(litho_strong_zones.x_right_strong.value)", level=3)
    end
end

function conditions_for_update(model::ModelData)::Bool
    iuse_strong_zones = model.geometry.parameters.litho_strong_zones.iuse_strong_zones.value
    if check_for_ttype_refinement(model) && ttype_refinement_parameters_are_set(model) && iuse_strong_zones == 1
        return true
    end
    return false
end

function ttype_refinement_parameters_are_set(model::ModelData)::Bool
    xo_highres = model.grids.parameters.refinement.xo_highres.value
    xf_highres = model.grids.parameters.refinement.xf_highres.value
    if !isnan(xo_highres) && !isnan(xf_highres)
        return true
    end
    return false
end

end 