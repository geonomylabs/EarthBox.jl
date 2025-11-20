""" Key inputs for different cases.

Case : Vext : dT   : Ext_Min : Ext_Max : Damage? : Ddist   : Dfac   : Lflow_sub : dt_erupt : tol   : max_iter : max_displ : C_trans,sa : R_pelagic  : Node
------------------------------------------------------------------------------------------------------------------------------------------------------------
   0 :  0.2 :  100 :    0.06 :    0.50 :       1 :  2500.0 :   10.0 :     20000 :    75000 : 1e-04 :      100 :      0.30 :    1.0e-04 :        0.0 : 
   1 :  0.2 :  100 :    0.06 :    0.10 :       1 :  2500.0 :   10.0 :     20000 :    75000 : 1e-04 :      100 :      0.30 :    1.0e-04 :        0.0 : 
   2 :  0.2 :  100 :    0.06 :    0.30 :       1 :  2500.0 :   10.0 :     20000 :    75000 : 1e-04 :      100 :      0.30 :    1.0e-04 :        0.0 : 
   3 :  0.2 :  100 :    0.06 :    0.50 :       1 :  2500.0 :   10.0 :     20000 :    75000 : 1e-04 :      100 :      0.30 :    1.0e-03 :        0.0 : 
   4 :  0.2 :  100 :    0.06 :    0.50 :       1 :  2500.0 :   10.0 :     20000 :    75000 : 1e-04 :      100 :      0.30 :    1.0e-02 :        0.0 : 
   5 :  0.2 :  100 :    0.06 :    0.50 :       1 :  2500.0 :   10.0 :     20000 :    75000 : 1e-04 :      100 :      0.30 :    1.0e-04 :       0.05 : 
   6 :  0.2 :  100 :    0.06 :    0.50 :       1 :  2500.0 :   10.0 :     20000 :    75000 : 1e-04 :      100 :      0.30 :    1.0e-04 :        0.1 : 
   7 :  0.2 :  100 :    0.06 :    0.50 :       1 :  2500.0 :   10.0 :     20000 :    75000 : 1e-04 :      100 :      0.30 :    1.0e-04 :        0.2 : 

"""
from earthbox.input_tools.case_inputs_tools import CaseInputsKeys
from earthbox.input_tools.case_inputs_tools import CaseCollectionType
from earthbox.input_tools.case_inputs_tools import initialize_cases
from earthbox.input_tools import case_builders
from earthbox.input_tools import case_printers
from earthbox.rock_properties import RhoCpOptionNames
from earthbox.melt_model import SolidusModels
from earthbox.melt_model import LiquidusModels


def get_case_inputs(
        model_case_name: str
) -> CaseCollectionType:
    """ Get key inputs for model setup.

    Returns
    -------
    key_inputs_active: dict
        Dictionary of key inputs for different cases. Key inputs are those that 
        change with different cases.

    """
    solidus_names = SolidusModels().options.option_names
    liquidus_names = LiquidusModels().options.option_names
    rhocp_names = RhoCpOptionNames()

    keys = CaseInputsKeys()
    base_case = {
        keys.use_mumps: True,
        # Model Setup
        keys.model_duration_myr: 32.0, # Myr
        keys.output_time_step_years: 200_000.0, # years
        keys.thickness_sticky_meters: 10_000.0, # meters
        keys.model_height_meters: 160_000.0, # meters
        keys.model_width_meters: 500_000.0, # meters
        # T-type Refinement
        keys.xo_highres: 150_000.0, # meters
        keys.xf_highres: 350_000.0, # meters
        keys.yf_highres: 50_000.0, # meters
        keys.dx_highres: 200.0, # meters
        keys.dx_lowres: 2000.0, # meters
        keys.dy_highres: 200.0, # meters
        keys.dy_lowres: 2000.0, # meters
        keys.dx_marker: 25.0, # meters
        keys.dy_marker: 25.0, # meters
        keys.dx_topo: 200.0, # meters
        # Rock Properties
        keys.rhocp_model: rhocp_names.temperature_dependent_waples,
        keys.maximum_heat_capacity: 1200.0, # J/kg/K
        keys.latent_heat_mantle: 400_000.0, # J/kg
        keys.latent_heat_oceanic_crust: 400_000.0, # J/kg
        keys.mantle_solidus: solidus_names.peridotite_katz2003,
        keys.mantle_liquidus: liquidus_names.peridotite_katz2003,
        # Thicknesses
        keys.thickness_lithosphere_meters: 125_000.0, # meters
        keys.thickness_continental_crust_meters: 35_000.0, # meters
        keys.thickness_upper_continental_crust_meters: 22_000.0, # meters
        # Velocity Boundary Conditions
        keys.full_velocity_extension_cm_yr: 0.2, # cm/yr
        keys.full_velocity_extension_cm_yr_step1: 0.6, # cm/yr
        keys.full_velocity_extension_cm_yr_step2: 2.0, # cm/yr
        keys.velocity_step1_time_myr: 8.0, # Myr
        keys.velocity_step2_time_myr: 16.0, # Myr
        # Temperature Boundary Conditions and Initial Temperature
        # T_bottom = 1345.0 + 0.4*25.0 = 1355.0
        # T_potential = 1355.0 - 0.4*150.0 = 1295.0
        keys.temperature_base_lithosphere_celsius: 1345.0, # degrees Celsius
        keys.use_bottom_transient: True,
        keys.delta_temperature_transient: 100.0, # degrees Celsius
        keys.start_time_bottom_transient_myr: 1.0, # Myr
        keys.end_time_bottom_transient_myr: 40.0, # Myr
        keys.adiabatic_temperature_gradient: 0.4, # C/km
        keys.radiogenic_heat_production_felsic_crust: 1.8e-6, # W/m^3
        keys.radiogenic_heat_production_mafic_crust: 0.5e-6, # W/m^3,
        keys.amplitude_perturbation_meters: 0.0, # meters
        # Stokes-continuity solver
        keys.maximum_number_of_global_plasticity_iterations: 100,
        keys.tolerance_picard: 1e-4,
        # Advection
        keys.marker_cell_displ_max: 0.3,
        # Hydrothermal Circulation
        keys.hydrothermal_nusselt_number: 2.0,
        keys.hydrothermal_max_depth_meters: 4000.0, # meters
        keys.use_plastic_strain_rate_for_hydrothermal: True,
        keys.hydrothermal_plastic_strain_rate_reference: 1e-14, # 1/s
        keys.use_plastic_strain_for_hydrothermal: True,
        keys.hydrothermal_plastic_strain_reference: 0.5,
        keys.hydrothermal_decay_length: 25_000.0, # meters
        keys.hydrothermal_buffer_distance: 25_000.0, # meters
        keys.sediment_thickness_threshold: 2500.0, # meters
        # Sediment Transport Model
        keys.precipitation_rate: 1.0, # m/yr
        keys.subaerial_transport_coefficient: 1.0e-4, # dimensionless
        keys.subaerial_slope_diffusivity: 0.25, # m2/yr
        keys.submarine_slope_diffusivity: 100.0, # m2/yr
        keys.submarine_diffusion_decay_depth: 1000.0, # meters
        keys.pelagic_sedimentation_rate: 0.0, # mm/yr
        # MeltDamageModel
        keys.use_melt_damage: True,
        keys.melt_damage_distance: 2_500.0, # meters
        keys.melt_damage_factor: 10.0,
        keys.use_probabilistic_melt_damage: True,
        keys.maximum_damage_probability: 0.8,
        keys.intermediate_damage_probability: 0.1,
        keys.magmatic_crust_height_threshold: 500.0, # meters
        keys.magmatic_crust_height_minimum: 750.0, # meters
        keys.magmatic_crust_height_intermediate: 2_000.0, # meters
        keys.magmatic_crust_height_maximum: 3_000.0, # meters
        # Lava Flow Model
        keys.characteristic_flow_length_submarine: 2_000.0, # meters
        keys.characteristic_flow_length_subaerial: 20_000.0, # meters
        keys.residual_lava_thickness_submarine: 30.0, # meters
        keys.residual_lava_thickness_subaerial: 30.0, # meters
        keys.porosity_initial_lava_flow: 0.0,
        # Extrusion Model
        keys.extrusion_volume_factor: 0.06,
        keys.extrusion_volume_factor_max: 0.5,
        keys.characteristic_magmatic_crust_height_min: 6_000.0, # meters
        keys.characteristic_magmatic_crust_height_max: 7_500.0, # meters
        keys.width_eruption_domain_fixed: 2_500.0, # meters, was 15_000.0
        keys.width_eruption_domain_fixed_max: 2_500.0, # meters, was 15_000.0
        keys.use_normal_eruption_location: True,
        keys.use_random_eruption_location: False,
        keys.use_eruption_interval: True,
        keys.eruption_interval_yr: 75_000.0, # years
        # Base Level Model
        keys.base_level_shift_meters: 5000.0, # meters
        keys.base_level_shift_end_time_myr: 16.0, # Myr
        # Rheological Parameters
        keys.plastic_healing_rate: 0.0, # 1/s
        # Serpentinization Model
        keys.maximum_serpentinization_depth_meters: 4_000.0, # meters,
        keys.maximum_serpentinization_rate: 1e-11 # 1/s
        }
    case_inputs = initialize_cases(base_case)

    # Variable maximum extrusion volume factor
    case_id = case_builders.build_case_group(
        case_inputs,
        case_id_ini=1,
        target_key=keys.extrusion_volume_factor_max,
        values = [0.1, 0.3]
        )
    # Variable subaerial transport coefficient
    case_id = case_builders.build_case_group(
        case_inputs,
        case_id_ini=case_id + 1,
        target_key=keys.subaerial_transport_coefficient,
        values = [1e-3, 1e-2],
        )
    # Variable pelagic sedimentation rate
    case_id = case_builders.build_case_group(
        case_inputs,
        case_id_ini=case_id + 1,
        target_key=keys.pelagic_sedimentation_rate,
        values = [0.05, 0.1, 0.2]
        )

    case_inputs_active = case_inputs[model_case_name]
    case_printers.print_case_info_sdr_formation(case_inputs, case_id)
    return case_inputs_active


def run_test():
    """ Run test
    """
    get_case_inputs('case0')


if __name__ == '__main__':
    run_test()
