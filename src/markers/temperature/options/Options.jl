module Options

import DataStructures: OrderedDict
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :TemperatureWave => 0,
    :Uniform => 1,
    :BoxConvection => 2,
    :Linear => 3,
    :AnalyticalThreeLayer => 4,
    :FourLinearSegments => 5,
    :FractureZone => 6,
    :HotBox => 7,
    :HalfSpaceSpreading => 8,
    :HalfSpaceDoubleSpreading => 9
)

const PDATA = get_eb_parameters()

const option_names = make_option_names(option_ids)

function get_options()::OrderedDict{Int, OptionState}
    return OrderedDict{Int, OptionState}(
        option_ids[option_names.TemperatureWave] =>
            OptionState(
                option_name=string(option_names.TemperatureWave),
                description="Simple temperature wave model whereby a square domain is rotated in a "
                *"circular velocity field. This model is used for the temperature wave benchmark.",
                required_parameters=[
                    "`$(PDATA.temperature_top.name)::Float64``: $(PDATA.temperature_top.description)",
                    "`$(PDATA.temperature_bottom.name)::Float64`: $(PDATA.temperature_bottom.description)",
                    "`$(PDATA.temperature_of_wave.name)::Float64`: $(PDATA.temperature_of_wave.description)",
                ],
                required_geometries=["None"]
            ),
        option_ids[option_names.Uniform] =>
            OptionState(
                option_name=string(option_names.Uniform),
                description="Uniform initial temperature temperature.",
                required_parameters=[
                    "`$(PDATA.temperature_uniform.name)::Float64`: $(PDATA.temperature_uniform.description)",
                ],
                required_geometries=["None"]
            ),
        option_ids[option_names.BoxConvection] =>
            OptionState(
                option_name=string(option_names.BoxConvection),
                description="Initial temperature model for convection in a box benchmark",
                required_parameters=[
                    "`$(PDATA.temperature_top.name)::Float64`: $(PDATA.temperature_top.description)",
                    "`$(PDATA.temperature_bottom.name)::Float64`: $(PDATA.temperature_bottom.description)",
                ],
                required_geometries=["None"]
            ),
        option_ids[option_names.Linear] =>
            OptionState(
                option_name=string(option_names.Linear),
                description="Initial temperature model using a linear increase with depth.",
                required_parameters=[
                    "`$(PDATA.temperature_top.name)::Float64`: $(PDATA.temperature_top.description)",
                    "`$(PDATA.temperature_bottom.name)::Float64`: $(PDATA.temperature_bottom.description)",
                ],
                required_geometries=["None"]
            ),
        option_ids[option_names.AnalyticalThreeLayer] =>
            OptionState(
                option_name=string(option_names.AnalyticalThreeLayer),
                description="Analytical steady-state three layer thermal model used to initialize marker temperature. " 
                *"The three layers are (1) the upper crust, (2) lower crust and (3) mantle "
                *"lithosphere all of which have user defined thermal conductivities and "
                *"radiogenic heat production. Thickness of the sticky layer above the three layers "
                *"and the thicknesses of the upper crust and lower crust are defined "
                *"during the initialization of earth layering and sticky geometry. Any "
                *"mantle lithosphere thicknesses defined during initialization are ignored "
                *"and the thickness of the thermal lithosphere (`thick_thermal_lithosphere`) "
                *"is used instead. The temperature at the top of the crust is set equal to the "
                *"temperature boundary condition applied at the top of the model (`temperature_top`) "
                *"The temperature at the base of the thermal lithosphere is set equal to the "
                *"`temperature_base_lith` parameter. If these parameters were set during boundary "
                *"condition initialization they do not need to be set here. Temperature below the "
                *"lithosphere is defined using an adiabatic gradient (`adiabatic_gradient`). A "
                *"thermal perturbation is also included in the model that involves a triangular "
                *"perturbation to the thickness of the thermal lithosphere. Note that if these "
                *"parameters are set they will be used to define the geotherm for the reference "
                *"lithosphere during reference lithosphere initialization.",
                required_parameters=[
                    "`$(PDATA.temperature_top.name)::Float64` : $(PDATA.temperature_top.description)",
                    "`$(PDATA.temperature_base_lith.name)::Float64` : $(PDATA.temperature_base_lith.description)",
                    "`$(PDATA.conductivity_upper_crust.name)::Float64` : $(PDATA.conductivity_upper_crust.description)",
                    "`$(PDATA.conductivity_lower_crust.name)::Float64` : $(PDATA.conductivity_lower_crust.description)",
                    "`$(PDATA.conductivity_mantle.name)::Float64` : $(PDATA.conductivity_mantle.description)",
                    "`$(PDATA.heat_production_upper_crust.name)::Float64` : $(PDATA.heat_production_upper_crust.description)",
                    "`$(PDATA.heat_production_lower_crust.name)::Float64` : $(PDATA.heat_production_lower_crust.description)",
                    "`$(PDATA.heat_production_mantle.name)::Float64` : $(PDATA.heat_production_mantle.description)",
                    "`$(PDATA.thick_thermal_lithosphere.name)::Float64` : $(PDATA.thick_thermal_lithosphere.description)",
                    "`$(PDATA.adiabatic_gradient.name)::Float64` : $(PDATA.adiabatic_gradient.description)",
                    "`$(PDATA.amplitude_perturbation.name)::Float64` : $(PDATA.amplitude_perturbation.description)",
                    "`$(PDATA.width_perturbation.name)::Float64` : $(PDATA.width_perturbation.description)",
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                ]

            ),
        option_ids[option_names.FourLinearSegments] =>
            OptionState(
                option_name=string(option_names.FourLinearSegments),
                description="Simple lithosphere geotherm composed of four linear segments with an "
                *"optional thermal perturbation located in the middle of the model. The "
                *"four segment are (1) sticky air/water, (2) continental crust, (3) mantle "
                *"lithosphere and (4) asthenosphere. The layer thicknesses are defined during "
                *"the initialization of earth layering and sticky air geometry parameters. "
                *"Sub-layering in the crust and mantle is lumped into crustal and mantle "
                *"lithosphere layers. Temperature in the sticky domain is constant and equal "
                *"to the temperature at the top of the model defined during boundary condition "
                *"initialization. Temperature within each additional layer is calculated "
                *"assuming a linear gradient defined by the temperatures set using this "
                *"input dictionary. A thermal perturbation is also included in the model "
                *"that involves a triangular perturbation to the depth of the asthenosphere.",
                required_parameters=[
                    "`$(PDATA.temperature_surface.name)::Float64` : $(PDATA.temperature_surface.description)",
                    "`$(PDATA.temperature_moho.name)::Float64` : $(PDATA.temperature_moho.description)",
                    "`$(PDATA.temperature_base_lith.name)::Float64` : $(PDATA.temperature_base_lith.description)",
                    "`$(PDATA.temperature_bottom.name)::Float64` : $(PDATA.temperature_bottom.description)",
                    "`$(PDATA.amplitude_perturbation.name)::Float64` : $(PDATA.amplitude_perturbation.description)",
                    "`$(PDATA.width_perturbation.name)::Float64` : $(PDATA.width_perturbation.description)",
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                ]
            ),
        option_ids[option_names.FractureZone] =>
            OptionState(
                option_name=string(option_names.FractureZone),
                description="Initial temperature condition for fracture zone model.",
                required_parameters=[
                    "`$(PDATA.temperature_top.name)::Float64` : $(PDATA.temperature_top.description)",
                    "`$(PDATA.temperature_bottom.name)::Float64` : $(PDATA.temperature_bottom.description)",
                    "`$(PDATA.age_lithosphere_left.name)::Float64` : $(PDATA.age_lithosphere_left.description)",
                    "`$(PDATA.age_lithosphere_right.name)::Float64` : $(PDATA.age_lithosphere_right.description)",
                    "`$(PDATA.adiabatic_gradient.name)::Float64` : $(PDATA.adiabatic_gradient.description)",
                    "`$(PDATA.thermal_lithosphere_depth_left.name)::Float64` : $(PDATA.thermal_lithosphere_depth_left.description)",
                    "`$(PDATA.thermal_lithosphere_depth_right.name)::Float64` : $(PDATA.thermal_lithosphere_depth_right.description)",
                    "`$(PDATA.thermal_diffusivity.name)::Float64` : $(PDATA.thermal_diffusivity.description)",
                ],
                required_geometries = [
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                    "[`FractureZone.initialize!`](@ref EarthBox.MaterialGeometry.FractureZone.initialize!)"
                ]
            ),
        option_ids[option_names.HotBox] =>
            OptionState(
                option_name=string(option_names.HotBox),
                description="Initial temperature condition for a hot box in a cold medium.",
                required_parameters=[
                    "`$(PDATA.temperature_top.name)::Float64` : $(PDATA.temperature_top.description)",
                    "`$(PDATA.temperature_bottom.name)::Float64` : $(PDATA.temperature_bottom.description)",
                ],
                required_geometries = [
                    "Geometric parameters are hard coded in the model. See `HotBox.jl` for more information."
                ]
            ),
        option_ids[option_names.HalfSpaceSpreading] =>
            OptionState(
                option_name=string(option_names.HalfSpaceSpreading),
                description="Initial temperature condition for half-space spreading model.",
                required_parameters=[
                    "`$(PDATA.temperature_top.name)::Float64` : $(PDATA.temperature_top.description)",
                    "`$(PDATA.adiabatic_gradient.name)::Float64` : $(PDATA.adiabatic_gradient.description)"
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                ]
            ),
        option_ids[option_names.HalfSpaceDoubleSpreading] =>
            OptionState(
                option_name=string(option_names.HalfSpaceDoubleSpreading),
                description="Initial temperature condition for half-space double spreading model.",
                required_parameters=[
                    "`$(PDATA.adiabatic_gradient.name)::Float64` : $(PDATA.adiabatic_gradient.description)"
                ],
                required_geometries = [
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                ]

            )
    )
end

end # module 