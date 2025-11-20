module Options

import DataStructures: OrderedDict
import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :solid_body_rotation => 0,
    :channel_flow_variable_conductivity => 1,
    :couette_flow_viscous_heating => 2,
    :channel_flow_non_steady_temperature => 3,
    :channel_flow_non_newtonian => 4,
    :rayleigh_taylor_instability => 5,
    :elastic_slab => 6,
    :viscoelastic_stress_buildup => 7,
    :box_convection_isoviscous_1a => 8,
    :plasticity_benchmark_kaus10 => 9,
    :seafloor_spreading => 10,
    :flexure_triangular_hole => 11,
    :viscoelastic_extension => 12,
    :viscoelastic_extension_asymmetric => 13,
    :viscoelastic_extension_depth_dependent => 14,
    :viscoelastic_extension_inflow_and_outflow_along_sides => 15,
    :viscoelastic_contraction => 16,
    :viscoelastic_contraction_asymmetric => 17,
    :simple_sedimentation => 18
)

const option_names = make_option_names(option_ids)

function get_options()::OrderedDict{Int, OptionState}
    OrderedDict{Int, OptionState}(
        option_ids[option_names.solid_body_rotation] =>
            OptionState(
                option_name=string(option_names.solid_body_rotation),
                description="Solid body rotation of a 10 km wide square region that " *
                "is 500 K hotter than the surrounding medium. Maximum velocity is " *
                "3.1536 cm/yr at a radius of 25 km centered in a 50 km x 50 km box. " *
                "Plots show results after 1 rotation (blue circles) and after 3 " *
                "rotation (black crosses) for cases with without subgrid diffusion " *
                "and thermal conductivity set to 1e-6 W/m/K. This benchmark tests " *
                "the code for numerical diffusion associated with the marker-in-cell " *
                "advection scheme. Resolution is 51 x 51 nodes with 25 randomly " *
                "distributed markers per cell. See Gerya, (2010), pg. 253.",
                time_steps=[1, 91, 181, 271],
                y_index=11
            ),
        option_ids[option_names.channel_flow_variable_conductivity] =>
            OptionState(
                option_name=string(option_names.channel_flow_variable_conductivity),
                description="Vertical constant-viscosity channel flow with temperature- " *
                "dependent thermal conductivity and shear heating. Boundary conditions are defined "*
                "using the `:VerticalChannelFlowShearHeating` boundary condition model type (See "*
                "[EarthBox.BoundaryConditions.initialize!](@ref)). "*
                "This benchmark tests the codes ability to handle strong thermal "*
                "conductivity gradients and shear heating terms in the heat equation. "*
                "Model parameters are as follows: ``L`` = 30 km, ``H`` = 11.25 km, ``eta`` = ``10^{19}`` Pa.s, "*
                "``P_beg`` = 3e7 Pa, ``P_end`` = 0, ``T_o`` = 298 K, ``k_o`` = 8 W/m/K, ``b`` = 1. "*
                "Resolution is 51 x 11 nodes with 25 randomly distributed markers per cell. "*
                "For details on the analytical solution see [gerya2010](@citet), pg. 254-255.",
                time_steps=[1, 11, 21, 31, 41, 51, 61, 71, 81, 91, 100],
                y_index=11
            ),
        option_ids[option_names.couette_flow_viscous_heating] =>
            OptionState(
                option_name=string(option_names.couette_flow_viscous_heating),
                description="Vertical Couette flow (simple shear in laterally limited " *
                "planer zone) with temperature-dependent Newtonian rheology and viscous " *
                "shear heating. Boundary conditions are defined using the `:VerticalCouetteFlow` "*
                "boundary condition model type (See [EarthBox.BoundaryConditions.initialize!](@ref)). "*
                "This benchmark tests how the code handles thermomechanical " *
                "effects of shear heating with temperature dependent viscosity. Model "*
                "parameters are as follows: ``L`` = 30 km, ``H`` = 11.25 km, ``A`` = ``10^{15}`` Pa.s, "*
                "``E_a`` = 150 kJ/mol, ``k`` = 2 W/m/K, ``T_o`` = 1000 K. Resolution is 51 x 11 nodes "*
                "with 25 randomly distributed markers per cell. For details on the analytical solution "*
                "see [gerya2010](@citet), pg. 250-253.",
                time_steps=[100],
                y_index=11
            ),
        option_ids[option_names.channel_flow_non_steady_temperature] =>
            OptionState(
                option_name=string(option_names.channel_flow_non_steady_temperature),
                description="Vertical Newtonian constant-viscosity channel flow with " *
                "non-steady temperature distribution. Boundary conditions are defined " *
                "using the `:VerticalChannelFlowIsoviscous` boundary condition " *
                "model type (See [EarthBox.BoundaryConditions.initialize!](@ref)). " *
                "The initial condition is defined by a constant vertical temperature " *
                "gradient and a constant temperature along bottom boundary. This benchmark tests the " *
                "codes ability to handle coupled advective and conductive terms in the " *
                "heat equation. Model parameters are as follows: ``L`` = 30 km, ``H`` = 11.25 km, " *
                "``eta`` = ``10^{19}`` Pa.s, ``P_{beg}`` = 1e5 Pa, ``P_{end}`` = 0, ``dT/dy`` = 40 K/km, "*
                "``T(y=0)`` = 1000 K. Resolution is 51 x 11 nodes with 25 randomly distributed markers per cell. " *
                "For details on the analytical solution see [gerya2010](@citet), pg. 247-250. "*
                "Note that a constant in equation 16.9 in [gerya2010](@citep) had to be modified from 8 " *
                "to 4 in order to match figure 16.5 in [gerya2010](@citep).",
                time_steps=[1, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 140, 160, 180, 200, 220, 240],
                y_index=11
            ),
        option_ids[option_names.channel_flow_non_newtonian] =>
            OptionState(
                option_name=string(option_names.channel_flow_non_newtonian),
                description="Vertical channel flow with non-Newtonian rheology. " *
                "Boundary conditions are defined using the `:VerticalChannelFlowNonNewtonian` boundary condition " *
                "model type (See [EarthBox.BoundaryConditions.initialize!](@ref)). " *
                "This benchmark tests the codes ability to handle flows where effective viscosity is " *
                "dependent on a power-law relationship with stress/strain rate. Model parameters are as follows: " *
                "``L`` = 10 km, ``H`` = 9.5 km, ``n`` = 3, ``C1`` = ``10^{-37}`` Pa^{-3}s^{-1}, ``P_{beg}`` = ``10^9`` Pa, " *
                "and ``P_{end}`` = 0. Resolution is 51 x 21 nodes with 25 randomly distributed markers per cell. " *
                "markers per cell. For details on the analytical solution see [gerya2010](@citet), pg. 246-247.",
                time_steps=[60],
                y_index=11
            ),
        option_ids[option_names.rayleigh_taylor_instability] =>
            OptionState(
                option_name=string(option_names.rayleigh_taylor_instability),
                description="Two-layer periodic Rayleigh-Taylor instability with sharp " *
                "contrast between density and viscosity fields. Constant viscosity and " *
                "density is defined with each layer. Boundary conditions are no-slip " *
                "along horizontal boundaries and symmetry along vertical walls. This " *
                "benchmark test the codes ability to calculate velocity for gravity " *
                "driven flows. Resolution is 51 x 51 nodes with 25 randomly distributed " *
                "markers per cell. See [gerya2010](@citet), pg. 242-244.",
                time_steps=[2],
                y_index=-1
            ),
        option_ids[option_names.elastic_slab] =>
            OptionState(
                option_name=string(option_names.elastic_slab),
                description="Recovery of an elastic slab that is surrounded by a " *
                "low-density, lower viscosity and higher shear modulus medium and " *
                "attached to the left vertical boundary. The slab is deformed under a " *
                "vertical gravity field (10 m/s/s) that is switched off after 20 Kyr to " *
                "allow recovery. The visco-elastic slab has a density of 4000 kg/m3, " *
                "viscosity of ``10^{27}`` Pa.s and a shear modulus of ``10^{10}`` Pa. The surrounding " *
                "visco-elastic medium has a density of 1 kg/m3, viscosity of ``10^{21}`` Pa.s, " *
                "shear modulus of ``10^{20}`` Pa. Boundary conditions are no-slip along the " *
                "left vertical boundary and free slip along all other boundaries. This " *
                "benchmark tests the codes ability to handle visco-elastic deformation " *
                "and conservation of elastic stress. Resolution is 51 x 51 nodes with 25 " *
                "randomly distributed markers per cell. For details on the analytical " *
                "solution see [gerya2010](@citet), pg. 261-263.",
                time_steps=[80, 100, 120, 140, 160],
                y_index=-1
            ),
        option_ids[option_names.viscoelastic_stress_buildup] =>
            OptionState(
                option_name=string(option_names.viscoelastic_stress_buildup),
                description="Stress buildup in an incompressible visco-elastic Maxwell " *
                "body undergoing uniform pure shear with constant strain rate. Boundary " *
                "conditions involve constant inward velocity along the horizontal " *
                "boundaries and constant outward velocity along the vertical boundaries. " *
                "Velocity boundary conditions are defined as a function of a constant " *
                "strain rate. This benchmark tests the codes ability to model elastic " *
                "buildup and the transition from primarily elastic deformation to " *
                "primarily viscous deformation. Model parameters: viscosity is 1e21 Pas, " *
                "shear modulus is 1e10 Pa, strain rate is 1e-14s-1. Resolution is 51 x 51 " *
                "nodes with 16 randomly distributed markers per cell. See Gerya, (2010), " *
                "pg. 260-261.",
                time_steps=[25],
                y_index=-1
            ),
        option_ids[option_names.box_convection_isoviscous_1a] =>
            OptionState(
                option_name=string(option_names.box_convection_isoviscous_1a),
                description="Stationary convection within constant viscosity aspect ratio " *
                "1 box. Free slip boundary conditions are applied along all boundaries. " *
                "Temperature boundary conditions are applied along the top and bottom " *
                "boundaries. Reflecting temperature boundary conditions are applied along " *
                "the vertical boundaries. This benchmark test the codes ability to handle " *
                "coupled thermal diffusion and convection for gravity driven flows. " *
                "Resolution is 51 x 51 nodes with 25 randomly distributed markers per " *
                "cell and local refinement along boundaries. See Gerya (2010), " *
                "pg. 255-260.",
                time_steps=[100],
                y_index=-1
            ),
        option_ids[option_names.plasticity_benchmark_kaus10] =>
            OptionState(
                option_name=string(option_names.plasticity_benchmark_kaus10),
                description="Extension benchmark similar to Kaus (2010) with visco-elasto " *
                "plastic rheology, single weak seed and free slip boundary condition " *
                "along the base of the model.",
                time_steps=[10],
                y_index=-1
            ),
        option_ids[option_names.flexure_triangular_hole] =>
            OptionState(
                option_name=string(option_names.flexure_triangular_hole),
                description="Flexure in response to a hole in a viscoelastic crustal layer.",
                time_steps=[60],
                y_index=-1
            ),
        option_ids[option_names.seafloor_spreading] =>
            OptionState(
                option_name=string(option_names.seafloor_spreading),
                description="Seafloor spreading via extension of lithospheric plate " *
                "benchmark with melting and extraction at ridge.",
                time_steps=[60],
                y_index=-1
            ),
        option_ids[option_names.viscoelastic_extension] =>
            OptionState(
                option_name=string(option_names.viscoelastic_extension),
                description="Simple extension of a viscoelastic layer undergoing water- " *
                "loaded symmetric extension with inflow along the top and bottom " *
                "boundaries and constant outflow along side boundaries. Velocity " *
                "inflow/outflow boundary conditions turned off to allow hydrostatic " *
                "equilibration. The numerical solution is compared to an analytical " *
                "solution for water-loaded local isostatic bathymetry.",
                time_steps=[1, 5, 10, 15, 20, 25, 30, 35, 40],
                y_index=-1
            ),
        option_ids[option_names.viscoelastic_extension_asymmetric] =>
            OptionState(
                option_name=string(option_names.viscoelastic_extension_asymmetric),
                description="Simple extension of a viscoelastic layer undergoing water- " *
                "loaded asymmetric extension with inflow along the top and bottom " *
                "boundaries. Velocity inflow/outflow boundary conditions turned off to " *
                "allow hydrostatic equilibration. The numerical solution is compared to " *
                "an analytical solution for water-loaded local isostatic bathymetry.",
                time_steps=[1, 5, 10, 15, 20, 25, 30, 35, 40],
                y_index=-1
            ),
        option_ids[option_names.viscoelastic_extension_depth_dependent] =>
            OptionState(
                option_name=string(option_names.viscoelastic_extension_depth_dependent),
                description="Simple extension of a viscoelastic layer undergoing water- " *
                "loaded symmetric depth-dependent extension with inflow along the top " *
                "and bottom boundaries. Velocity inflow/outflow boundary conditions turned " *
                "off to allow hydrostatic equilibration. The numerical solution is " *
                "compared to an analytical solution for water-loaded local isostatic " *
                "bathymetry.",
                time_steps=[1, 5, 10, 15, 20, 25, 30, 35, 40],
                y_index=-1
            ),
        option_ids[option_names.viscoelastic_extension_inflow_and_outflow_along_sides] =>
            OptionState(
                option_name=string(option_names.viscoelastic_extension_inflow_and_outflow_along_sides),
                description="Simple extension of a viscoelastic layer undergoing water- " *
                "loaded symmetric extension with inflow along side and top boundaries " *
                "and outflow along side boundaries. Velocity inflow/outflow boundary " *
                "conditions along side boundaries are turned off to allow hydrostatic " *
                "equilibration. The numerical solution is compared to an analytical " *
                "solution for water-loaded local isostatic bathymetry.",
                time_steps=[1, 5, 10, 15, 20, 25, 30, 35, 40],
                y_index=-1
            ),
        option_ids[option_names.viscoelastic_contraction] =>
            OptionState(
                option_name=string(option_names.viscoelastic_contraction),
                description="Simple contraction of a viscoelastic layer undergoing water- " *
                "loaded symmetric contraction with outflow along the top and bottom " *
                "boundaries and constant inflow along side boundaries. Velocity " *
                "inflow/outflow boundary conditions are turned off to allow hydrostatic " *
                "equilibration. The numerical solution is compared to an analytical " *
                "solution for water-loaded local isostatic bathymetry.",
                time_steps=[1, 2, 3, 4, 5, 6, 7, 8],
                y_index=-1
            ),
        option_ids[option_names.viscoelastic_contraction_asymmetric] =>
            OptionState(
                option_name=string(option_names.viscoelastic_contraction_asymmetric),
                description="Simple contraction of a viscoelastic layer undergoing water- " *
                "loaded asymmetric contraction with outflow along the top and bottom " *
                "boundaries and constant inflow along side boundaries. Velocity " *
                "inflow/outflow boundary conditions are turned off to allow hydrostatic " *
                "equilibration. The numerical solution is compared to an analytical " *
                "solution for water-loaded local isostatic bathymetry.",
                time_steps=[1, 2, 3, 4, 5, 6, 7, 8],
                y_index=-1
            ),
        option_ids[option_names.simple_sedimentation] =>
            OptionState(
                option_name=string(option_names.simple_sedimentation),
                description=(
                    "Sedimentation with compaction on a mantle substrate with a constant pelagic "
                    *"sedimentation rate of 20 mm/yr. "
                    *"The model domain is 140 km x 140 km with 251 x 71 nodes and 16 markers per "
                    *"cell in the horizontal and vertical directions. "
                    *"The topography marker chain has 561 uniformly spaced nodes. "
                    *"The model is initialized with a 10 km sticky layer with 5 km of sticky water "
                    *"and a constant sea level is used during the simulation. "
                    *"The Stokes-continuity and heat equation are not solved. "
                    *"This benchmark tests the codes ability to model sedimentation and compaction "
                    *"against a known analytical solution. "
                ),
                time_steps=[5],
                y_index=-1
            )
    )
end

end # module 