module EarthLayering

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

export initialize!

const PDATA = get_eb_parameters()

struct ValidInputNames
    thick_upper_crust::Symbol
    thick_lower_crust::Symbol
    thick_upper_lith::Symbol
    thick_middle_lith::Symbol
    thick_lower_lith::Symbol
    thick_lith::Symbol
    thick_crust::Symbol
    thick_mantle_lith::Symbol
end

""" Initialize Earth layering geometry.

# Arguments
- model::ModelData: The model data container used to store model parameters and arrays.

# Keyword arguments:
- `thick_lith::Float64`: 
    - $(PDATA.thick_lith.description)
- `thick_crust::Float64`: 
    - $(PDATA.thick_crust.description)
- `thick_upper_crust::Float64`: 
    - $(PDATA.thick_upper_crust.description)
- `thick_lower_crust::Float64`: 
    - $(PDATA.thick_lower_crust.description)
- `thick_upper_lith::Float64`: 
    - $(PDATA.thick_upper_lith.description)
- `thick_middle_lith::Float64`: 
    - $(PDATA.thick_middle_lith.description)
- `thick_lower_lith::Float64`: 
    - $(PDATA.thick_lower_lith.description)

# Returns
- Nothing

# Initialization Scenarios

**Scenario 1**: Layer thicknesses provided as input are used as-is for the upper continental crust 
   (`thick_upper_crust`), lower continental crust (`thick_lower_crust`), upper mantle lithosphere 
   (`thick_upper_lith`), middle mantle lithosphere (`thick_middle_lith`) and lower mantle 
   lithosphere (`thick_lower_lith`).

**Scenario 2**: If the total thickness of the lithosphere (`thick_lith`), 
the total thickness of the continental crust (`thick_crust`), and the thickness
of the upper continental crust (`thick_upper_crust`) are provided as inputs, remaining
thicknesses are derived as follows:
- The thickness of the lower continental crust (`thick_lower_crust`) will be calculated as the 
   difference between the total thickness of the continental crust (`thick_crust`) and the thickness 
   of the upper continental crust (`thick_upper_crust`).
- The thickness of the mantle lithosphere (`thick_mantle_lith`) will be calculated as the difference 
   between the thickness of the total lithosphere (`thick_lith`) and the thickness of the continental 
   crust (`thick_crust`).
- The thickness of the upper mantle lithosphere (`thick_upper_lith`) will be set to one-third of the 
   thickness of the mantle lithosphere (`thick_mantle_lith`).
- The thickness of the middle mantle lithosphere (`thick_middle_lith`) will be set to one-third of the 
   thickness of the mantle lithosphere (`thick_mantle_lith`).
- The thickness of the lower mantle lithosphere (`thick_lower_lith`) will be set to one-third of the 
   thickness of the mantle lithosphere (`thick_mantle_lith`).

# Example 1: Direct input of layer thicknesses

```julia
eb = EarthBoxState(
    xnum = 200, ynum = 100, xsize = 500000.0, ysize = 120000.0,
    dx_marker = 100.0, dy_marker = 100.0
    )
model = eb.model_manager.model # ModelData container
EarthBox.MaterialGeometry.EarthLayering.initialize!(
    model, 
    thick_upper_crust = 22_000.00, 
    thick_lower_crust = 18_000.0, 
    thick_upper_lith  = 28_333.3, 
    thick_middle_lith = 28_333.3, 
    thick_lower_lith  = 22_000.0
)
```

# Example 2: Calculating layer thicknesses from total lithosphere and crust thicknesses

```julia
eb = EarthBoxState(
    xnum = 200, ynum = 100, xsize = 500000.0, ysize = 120000.0,
    dx_marker = 100.0, dy_marker = 100.0
)
model = eb.model_manager.model # ModelData container
EarthBox.MaterialGeometry.EarthLayering.initialize!(
    model, 
    thick_lith        = 100_000.0, 
    thick_crust       = 40_000.0, 
    thick_upper_crust = 22_000.0
)
```

"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    update_parameters!(model)
    update_thick_lith(model)
    return nothing
end

function update_parameters!(model::ModelData)
    earth_layering = model.geometry.parameters.earth_layering
    
    thick_lith        = earth_layering.thick_lith.value
    thick_crust       = earth_layering.thick_crust.value
    thick_upper_crust = earth_layering.thick_upper_crust.value
    
    if parameters_are_loaded_for_derived_inputs(thick_lith, thick_crust, thick_upper_crust)
        thick_mantle_lith = calculate_mantle_lithosphere_thickness(thick_lith, thick_crust)
        thick_lower_crust = calculate_lower_crust_thickness(thick_crust, thick_upper_crust)
        thick_upper_lith = thick_middle_lith = thick_lower_lith = thick_mantle_lith / 3.0
        
        earth_layering.thick_lower_crust.value = thick_lower_crust
        earth_layering.thick_upper_lith.value  = thick_upper_lith
        earth_layering.thick_middle_lith.value = thick_middle_lith
        earth_layering.thick_lower_lith.value  = thick_lower_lith
        earth_layering.thick_mantle_lith.value = thick_mantle_lith

        print_info("Calculated derived parameters for Earth Layering", level=1)
        print_info("Inputs:", level=2)
        print_info("thick_lith: $(thick_lith)", level=3)
        print_info("thick_crust: $(thick_crust)", level=3)
        print_info("thick_upper_crust: $(thick_upper_crust)", level=3)
        print_info("Outputs:", level=2)
        print_info("thick_mantle_lith: $(thick_mantle_lith)", level=3)
        print_info("thick_lower_crust: $(thick_lower_crust)", level=3)
        print_info("thick_upper_lith: $(thick_upper_lith)", level=3)
        print_info("thick_middle_lith: $(thick_middle_lith)", level=3)
        print_info("thick_mantle_lith: $(thick_mantle_lith)", level=3)
    end
end

function parameters_are_loaded_for_derived_inputs(
    thick_lith::Float64,
    thick_crust::Float64,
    thick_upper_crust::Float64
)::Bool
    if isnan(thick_lith) || isnan(thick_crust) || isnan(thick_upper_crust)
        return false
    end
    print_info("Earth Layering parameters are loaded for derived parameters", level=1)
    print_info("thick_lith: $(thick_lith)", level=2)
    print_info("thick_crust: $(thick_crust)", level=2)
    print_info("thick_upper_crust: $(thick_upper_crust)", level=2)
    return true
end

function update_thick_lith(model::ModelData)::Nothing
    earth_layering = model.geometry.parameters.earth_layering
    thick_upper_lith = earth_layering.thick_upper_lith.value
    thick_middle_lith = earth_layering.thick_middle_lith.value
    thick_lower_lith = earth_layering.thick_lower_lith.value

    # If all three mantle lithosphere layers are defined, then update thick_lith
    check = false
    if !isnan(thick_upper_lith) && !isnan(thick_middle_lith) && !isnan(thick_lower_lith)
        check = true
    end
    if check
        thick_lith = thick_upper_lith + thick_middle_lith + thick_lower_lith
        earth_layering.thick_lith.value = thick_lith
        print_info("Updated thick_lith: $(thick_lith)", level=1)
    end
    return nothing
end

function calculate_lower_crust_thickness(
    thick_crust::Float64,
    thick_upper_crust::Float64
)::Float64
    return thick_crust - thick_upper_crust
end

function calculate_mantle_lithosphere_thickness(
    thick_lith::Float64,
    thick_crust::Float64
)::Float64
    return thick_lith - thick_crust
end

end 