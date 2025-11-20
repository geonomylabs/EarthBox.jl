module Options

import DataStructures: OrderedDict
import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names
import ..Registry

const option_ids = Dict{Symbol, Int}(
    :UniformComposition => 0,
    :RayleighTaylorInstabilityBenchmark => 1,
    :ViscousBlockBenchmark => 2,
    :SandboxShortening => 3,
    :SandboxExtension => 4,
    :ElasticSlab => 5,
    :LithosphericExtensionWeakSeed => 6,
    :LithosphericExtensionNoSeed => 7,
    :PlasticityBenchmark => 8,
    :PlasticityBenchmarkWeakNotchKaus10 => 9,
    :FlexureTriangularHole => 10,
    :LithosphericExtensionWeakFault => 11,
    :LithosphericExtensionLateralStrongZones => 12,
    :LithosphericExtensionMantlePlume => 13,
    :SeafloorSpreading => 14,
    :SimplePlume => 15,
    :FractureZone => 16,
    :SimpleSediment => 17,
    :StickyOnMantle => 18
)

const option_names = make_option_names(option_ids)

function get_options()::OrderedDict{Int, OptionState}
    domains = Registry.MaterialDomainsRegistry()
    
    return OrderedDict{Int, OptionState}(
        option_ids[option_names.UniformComposition] =>
            OptionState(
                option_name=string(option_names.UniformComposition),
                description="Material model type with uniform composition.",
                required_domains=[domains.general_domain],
                required_geometries=["None"]
            ),
        option_ids[option_names.RayleighTaylorInstabilityBenchmark] =>
            OptionState(
                option_name=string(option_names.RayleighTaylorInstabilityBenchmark),
                description="Material model type for the Rayleigh-Taylor instability benchmark from Ramberg.",
                required_domains=[
                    domains.layer1_rayleigh_taylor, 
                    domains.layer2_rayleigh_taylor
                ],
                required_geometries=["[`RayleighTaylor.initialize!`](@ref EarthBox.MaterialGeometry.RayleighTaylor.initialize!)"]
            ),
        option_ids[option_names.ViscousBlockBenchmark] =>
            OptionState(
                option_name=string(option_names.ViscousBlockBenchmark),
                description="Material model type for viscous block benchmark.",
                required_domains=[domains.general_domain],
                required_geometries=["None"]
            ),
        option_ids[option_names.SandboxShortening] =>
            OptionState(
                option_name=string(option_names.SandboxShortening),
                description="Material model type for sandbox shortening experiments.",
                required_domains=[
                    domains.atmosphere, 
                    domains.sand_a, 
                    domains.sand_b, 
                    domains.microbeads,
                    domains.mobile_wall
                ],
                required_geometries=[
                    "[`Sandbox.initialize!`](@ref EarthBox.MaterialGeometry.Sandbox.initialize!)",
                    "[`MobileWall.initialize!`](@ref EarthBox.MaterialGeometry.MobileWall.initialize!)"
                ]
            ),
        option_ids[option_names.SandboxExtension] =>
            OptionState(
                option_name=string(option_names.SandboxExtension),
                description="Material model type for sandbox extension experiments.",
                required_domains=[
                    domains.atmosphere, domains.sand_a, domains.sand_b, 
                    domains.pdms_layer, domains.mobile_wall, domains.plate_extension
                ],
                required_geometries=[
                    "[`Sandbox.initialize!`](@ref EarthBox.MaterialGeometry.Sandbox.initialize!)",
                    "[`MobileWall.initialize!`](@ref EarthBox.MaterialGeometry.MobileWall.initialize!)"
                ]
            ),
        option_ids[option_names.ElasticSlab] =>
            OptionState(
                option_name=string(option_names.ElasticSlab),
                description="Material model type for elastic slab benchmark. This material involves "
                *"an elastic slab with a checkerboard structure composed of two material " 
                *"domains called `ElasticSlabA` and `ElasticSlabB` surrounded by two fluid domains called "
                *"`FluidAroundElasticSlabA` and `FluidAroundElasticSlabB`. Materials A and B have equivalent "
                *"properties but different colors to visualize the checkerboard structure",
                required_domains=[
                    domains.fluid_around_elastic_slab_a, 
                    domains.fluid_around_elastic_slab_b, 
                    domains.elastic_slab_a, 
                    domains.elastic_slab_b
                ],
                required_geometries=["None"]
            ),
        option_ids[option_names.LithosphericExtensionWeakSeed] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionWeakSeed),
                description="Material model type for lithospheric extension similar to Liao et al. (2014)",
                required_domains=[
                    domains.atmosphere, 
                    domains.ocean, 
                    domains.upper_continental_crust,
                    domains.lower_continental_crust,
                    domains.upper_mantle_lithosphere, 
                    domains.middle_mantle_lithosphere, 
                    domains.lower_mantle_lithosphere,
                    domains.asthenosphere, 
                    domains.weak_seed_mantle
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                    "[`WeakSeed.initialize!`](@ref EarthBox.MaterialGeometry.WeakSeed.initialize!)"
                ]
            ),
        option_ids[option_names.LithosphericExtensionNoSeed] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionNoSeed),
                description="Material model type for lithospheric extension model similar to Brune et al. (2014).",
                required_domains=[
                    domains.atmosphere, domains.ocean, 
                    domains.upper_continental_crust,
                    domains.lower_continental_crust,
                    domains.upper_mantle_lithosphere,
                    domains.middle_mantle_lithosphere, 
                    domains.lower_mantle_lithosphere,
                    domains.asthenosphere
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)"
                ]
            ),
        option_ids[option_names.PlasticityBenchmark] =>
            OptionState(
                option_name=string(option_names.PlasticityBenchmark),
                description="Material model type for plasticity benchmark.",
                required_domains=[
                    domains.atmosphere,
                    domains.upper_continental_crust,
                    domains.lower_continental_crust,
                    domains.upper_mantle_lithosphere, 
                    domains.middle_mantle_lithosphere, 
                    domains.lower_mantle_lithosphere, 
                    domains.asthenosphere, 
                    domains.weak_seed_mantle
                    ],
                required_geometries=["None"]
            ),
        option_ids[option_names.PlasticityBenchmarkWeakNotchKaus10] =>
            OptionState(
                option_name=string(option_names.PlasticityBenchmarkWeakNotchKaus10),
                description="Material model type for plasticity benchmark of Kaus (2010) with sticky air, crust and weak notch",
                required_domains=[
                    domains.atmosphere, 
                    domains.upper_continental_crust, 
                    domains.lower_continental_crust, 
                    domains.upper_mantle_lithosphere, 
                    domains.middle_mantle_lithosphere, 
                    domains.lower_mantle_lithosphere, 
                    domains.asthenosphere, 
                    domains.weak_seed_mantle
                    ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                    "[`WeakSeed.initialize!`](@ref EarthBox.MaterialGeometry.WeakSeed.initialize!)"
                ]
            ),
        option_ids[option_names.FlexureTriangularHole] =>
            OptionState(
                option_name=string(option_names.FlexureTriangularHole),
                description="Material model type for flexure with a triangular hole benchmark.",
                required_domains=[
                    domains.atmosphere, 
                    domains.ocean, 
                    domains.upper_continental_crust, 
                    domains.lower_continental_crust, 
                    domains.upper_mantle_lithosphere, 
                    domains.middle_mantle_lithosphere,
                    domains.lower_mantle_lithosphere, 
                    domains.asthenosphere
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                    "[`CrustalHole.initialize!`](@ref EarthBox.MaterialGeometry.CrustalHole.initialize!)"
                ]
            ),
        option_ids[option_names.LithosphericExtensionWeakFault] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionWeakFault),
                description="Material model type for lithospheric extension using a weak fault similar to Davis et al. (2017).",
                required_domains=[
                    domains.atmosphere, 
                    domains.ocean, 
                    domains.upper_continental_crust,
                    domains.lower_continental_crust,
                    domains.upper_mantle_lithosphere,
                    domains.middle_mantle_lithosphere, 
                    domains.lower_mantle_lithosphere,
                    domains.asthenosphere,
                    domains.weak_mantle_fault_zone,
                    domains.weak_crustal_fault_zone
                    ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                    "[`WeakFault.initialize!`](@ref EarthBox.MaterialGeometry.WeakFault.initialize!)"
                ]
            ),
        option_ids[option_names.LithosphericExtensionLateralStrongZones] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionLateralStrongZones),
                description="Material model type for lithospheric extension using lateral strong "
                *"zones similar to [naliboff17](@cite).",
                required_domains=
                required_domains=[
                    domains.atmosphere,
                    domains.ocean, 
                    domains.upper_continental_crust,
                    domains.upper_continental_crust_strong_zone,
                    domains.lower_continental_crust, 
                    domains.lower_continental_crust_strong_zone,
                    domains.upper_mantle_lithosphere, 
                    domains.middle_mantle_lithosphere, 
                    domains.lower_mantle_lithosphere, 
                    domains.lithospheric_mantle_strong_zone,
                    domains.asthenosphere
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                    "[`LithoStrongZones.initialize!`](@ref EarthBox.MaterialGeometry.LithoStrongZones.initialize!)"
                ]
            ),
        option_ids[option_names.LithosphericExtensionMantlePlume] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionMantlePlume),
                description="Material model type for lithospheric extension with mantle plume.",
                required_domains=[
                    domains.atmosphere, domains.ocean, 
                    domains.sedimentary_basin, 
                    domains.upper_continental_crust, 
                    domains.upper_continental_crust_strong_zone, 
                    domains.lower_continental_crust, 
                    domains.lower_continental_crust_strong_zone, 
                    domains.upper_mantle_lithosphere, 
                    domains.middle_mantle_lithosphere, 
                    domains.lower_mantle_lithosphere, 
                    domains.lithospheric_mantle_strong_zone, 
                    domains.asthenosphere, 
                    domains.weak_seed_mantle, 
                    domains.mantle_plume
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)",
                    "[`LithoStrongZones.initialize!`](@ref EarthBox.MaterialGeometry.LithoStrongZones.initialize!)",
                    "[`WeakSeed.initialize!`](@ref EarthBox.MaterialGeometry.WeakSeed.initialize!)",
                    "[`Plume.initialize!`](@ref EarthBox.MaterialGeometry.Plume.initialize!)"
                ]
            ),
        option_ids[option_names.SeafloorSpreading] =>
            OptionState(
                option_name=string(option_names.SeafloorSpreading),
                description="Material model type for seafloor spreading benchmark.",
                required_domains=[
                    domains.atmosphere, 
                    domains.ocean, 
                    domains.mantle_lithosphere, 
                    domains.asthenosphere
                ],
                required_geometries=[
                    "[`EarthLayering.initialize!`](@ref EarthBox.MaterialGeometry.EarthLayering.initialize!)",
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)"
                ]
            ),
        option_ids[option_names.SimplePlume] =>
            OptionState(
                option_name=string(option_names.SimplePlume),
                description="Material model for a simple plume head.",
                required_domains=[
                    domains.atmosphere, 
                    domains.ocean, 
                    domains.asthenosphere, 
                    domains.mantle_plume
                ],
                required_geometries=[
                    "[`Plume.initialize!`](@ref EarthBox.MaterialGeometry.Plume.initialize!)"
                ]
            ),
        option_ids[option_names.FractureZone] =>
            OptionState(
                option_name=string(option_names.FractureZone),
                description="Material model type for a fracture zone.",
                required_domains=[
                    domains.atmosphere, domains.ocean, 
                    domains.sedimentary_basin, 
                    domains.basaltic_oceanic_crust, 
                    domains.gabbroic_oceanic_crust, 
                    domains.mantle_lithosphere, 
                    domains.weak_mantle_lithosphere, 
                    domains.asthenosphere
                ],
                required_geometries=[
                    "[`FractureZone.initialize!`](@ref EarthBox.MaterialGeometry.FractureZone.initialize!)"
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)"
                ]
            ),
        option_ids[option_names.SimpleSediment] =>
            OptionState(
                option_name=string(option_names.SimpleSediment),
                description="Material model for simple sediment on mantle setup.",
                required_domains=[
                    domains.atmosphere, 
                    domains.sedimentary_basin, 
                    domains.asthenosphere
                ],
                required_geometries=[
                    "[`Plume.initialize!`](@ref EarthBox.MaterialGeometry.Plume.initialize!)"
                ]
            ),
        option_ids[option_names.StickyOnMantle] =>
            OptionState(
                option_name=string(option_names.StickyOnMantle),
                description="Material model for sticky air/water on mantle.",
                required_domains=[domains.ocean, domains.asthenosphere],
                required_geometries=[
                    "[`StickyAirGeometry.initialize!`](@ref EarthBox.MaterialGeometry.StickyAirGeometry.initialize!)"
                ]
            )
    )
end

end # module 