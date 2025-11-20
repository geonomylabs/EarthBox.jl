module Options

import DataStructures: OrderedDict
import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :LithosphericExtensionMovingBoundaries => 0,
    :LithosphericExtensionFixedBoundaries => 1,
    :LithosphericExtensionFixedBoundariesAsymmetric => 2,
    :LithosphericExtensionDepthDependent => 3,
    :LithosphericExtensionInflowAndOutflowAlongSides => 4,
    :ExtensionNoBottomInflow => 5,
    :PureShearDeformation => 6,
    :VerticalCouetteFlow => 7,
    :SolidBodyRotation => 8,
    :RayleighTaylorInstability => 9,
    :BoxConvection => 10,
    :ElasticSlab => 11,
    :ViscousBlock => 12,
    :SandboxShortening => 13,
    :SandboxExtension => 14,
    :VerticalChannelFlowIsoviscous => 15,
    :VerticalChannelFlowNonNewtonian => 16,
    :VerticalChannelFlowShearHeating => 17,
    :FractureZoneContractionFixedBoundaries => 18,
    :FractureZoneContractionFixedBoundariesAsymmetric => 19,
    :ViscoelasticContraction => 20,
    :ViscoelasticContractionAsymmetric => 21
)

const option_names = make_option_names(option_ids)

function get_options()::OrderedDict{Int, OptionState}
    return OrderedDict{Int, OptionState}(
        option_ids[option_names.LithosphericExtensionMovingBoundaries] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionMovingBoundaries),
                description="This boundary condition type is used for extension models with moving boundaries. " *
                          "A fixed free slip upper boundary is used with moving lower and side boundaries. " *
                          "Orthogonal velocity is prescribed along side and bottom boundaries. " *
                          "Regridding is applied with no marker recycling since boundaries are moving. " *
                          "Temperature is prescribed on the top and bottom, and zero heat flux is applied " *
                          "along side boundaries.",
                bools=Dict{Symbol, Bool}(:extending_side_boundaries => true)
            ),
        option_ids[option_names.LithosphericExtensionFixedBoundaries] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionFixedBoundaries),
                description="This boundary condition type is for extension models with fixed boundaries. " *
                          "Free slip with orthogonal velocity is prescribed along top, bottom and side " *
                          "boundaries to ensure volume conservation. Markers are recycled with special " *
                          "treatment in the sticky air region to ensure constant sticky air volume. " *
                          "Temperature is prescribed on the top and bottom, and zero heat flux is applied " *
                          "along side boundaries.",
                bools=Dict{Symbol, Bool}(:stepping_velocity_bc => true)
            ),
        option_ids[option_names.LithosphericExtensionFixedBoundariesAsymmetric] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionFixedBoundariesAsymmetric),
                description="This boundary condition type is for extension models with fixed boundaries and " *
                          "asymmetric boundary conditions. Free slip in the y-direction and zero lateral " *
                          "orthogonal velocity are applied on the left boundary, and uniform lateral full " *
                          "extension velocity with free slip in y-direction is applied along the right " *
                          "boundary. Free slip with orthogonal velocity is prescribed along top and bottom " *
                          "boundaries to ensure volume conservation. Markers are recycled with special " *
                          "treatment in the sticky air region to ensure constant sticky air volume. " *
                          "Temperature is prescribed on the top and bottom, and zero heat flux is applied " *
                          "along side boundaries.",
                bools=Dict{Symbol, Bool}(:stepping_velocity_bc => true)
            ),
        option_ids[option_names.LithosphericExtensionDepthDependent] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionDepthDependent),
                description="This boundary condition type is used for extension models with fixed boundaries " *
                          "and side boundaries with free slip and depth-dependent prescribed orthogonal " *
                          "velocity. Along top and bottom boundaries free slip is applied with prescribed " *
                          "orthogonal velocity calculated to ensure volume conservation. Particles are " *
                          "recycled (no regridding) with special treatment in sticky air region to ensure " *
                          "constant sticky air volume. Temp. boundary conditions on top and bottom, zero " *
                          "flux on sides.",
                bools=Dict{Symbol, Bool}(:stepping_velocity_bc => true, :depth_dependent_velocity => true)
            ),
        option_ids[option_names.LithosphericExtensionInflowAndOutflowAlongSides] =>
            OptionState(
                option_name=string(option_names.LithosphericExtensionInflowAndOutflowAlongSides),
                description="This boundary condition type is used for extension models with fixed boundaries " *
                          "and side boundaries with free slip and depth-dependent prescribed orthogonal " *
                          "velocity for inflow and outflow. Along top and bottom boundary free slip is " *
                          "applied with prescribed orthogonal velocity calculated to ensure volume " *
                          "conservation in the sticky layer. Particles are recycled (no regridding) with " *
                          "special treatment in sticky and side inflow regions to ensure. Temperature " *
                          "boundary conditions on top and bottom, zero flux on sides.",
                bools=Dict{Symbol, Bool}(:stepping_velocity_bc => true)
            ),
        option_ids[option_names.ExtensionNoBottomInflow] =>
            OptionState(
                option_name=string(option_names.ExtensionNoBottomInflow),
                description="This boundary condition type is used for the simple extension with constant " *
                          "outflow along the side boundary, conservative inflow along the top boundary and " *
                          "a closed free slip boundary along the bottom. Velocity is prescribed along side " *
                          "and top boundaries to ensure the conservation of volume.",
                bools=Dict{Symbol, Bool}(:stepping_velocity_bc => true)
            ),
        option_ids[option_names.PureShearDeformation] =>
            OptionState(
                option_name=string(option_names.PureShearDeformation),
                description="This boundary conditions type is used for a pure shear experiment. Free slip " *
                          "is applied along boundaries with orthogonal prescribed velocity. Temperature is " *
                          "prescribed on top, bottom and side boundaries. No marker recycling is used.",
            ),
        option_ids[option_names.VerticalCouetteFlow] =>
            OptionState(
                option_name=string(option_names.VerticalCouetteFlow),
                description="Boundary conditions for vertical Couette flow (simple shear). Open channel " *
                          "boundary conditions are applied along the top and bottom boundaries to model " *
                          "vertical flow in an infinite channel. This is accomplished by prescribing zero " *
                          "vertical velocity gradients along top and bottom boundaries and by prescribing " *
                          "a vertical pressure gradient, which for this case is zero. No-slip bcs on are " *
                          "applied on side boundaries with prescribed y-velocity along the right boundary " *
                          "to induce shear flow. Zero vertical temperature gradients are prescribed along " *
                          "top and bottom boundaries. Prescribed temperature is applied along left " *
                          "boundary. An insulating boundary condition (zero lateral heat flux) is " *
                          "implemented along the right boundary. Marker recycling takes into account " *
                          "vertical channel flow.",
            ),
        option_ids[option_names.SolidBodyRotation] =>
            OptionState(
                option_name=string(option_names.SolidBodyRotation),
                description="This boundary condition type is used for the solid body rotation benchmark " *
                          "with fixed boundaries (Stokes equations are not solved). Prescribed thermal " *
                          "gradients are applied along side boundaries.",
            ),
        option_ids[option_names.RayleighTaylorInstability] =>
            OptionState(
                option_name=string(option_names.RayleighTaylorInstability),
                description="This boundary condition type is used for the Rayleigh-Taylor instability " *
                          "benchmark from Ramberg (1968). No slip boundary conditions are applied on top " *
                          "and bottom boundaries and free slip boundary condition are applied on side " *
                          "boundaries.",
            ),
        option_ids[option_names.BoxConvection] =>
            OptionState(
                option_name=string(option_names.BoxConvection),
                description="This boundary conditions type is used for the convection in a box benchmark. " *
                          "Free slip boundary conditions are applied on all sides. Temperature is " *
                          "prescribed on top and bottom boundaries. Side boundaries have insulating " *
                          "thermal boundary conditions.",
            ),
        option_ids[option_names.ElasticSlab] =>
            OptionState(
                option_name=string(option_names.ElasticSlab),
                description="This boundary conditions type is used for the elastic slab experiment. No " *
                          "slip boundary conditions are applied on the left boundary. Free slip boundary " *
                          "conditions are applied on top, bottom and right boundaries.",
            ),
        option_ids[option_names.ViscousBlock] =>
            OptionState(
                option_name=string(option_names.ViscousBlock),
                description="This boundary conditions type is used for the viscous block benchmark. Free " *
                          "slip boundary conditions are applied along all boundaries.",
            ),
        option_ids[option_names.SandboxShortening] =>
            OptionState(
                option_name=string(option_names.SandboxShortening),
                description="This boundary condition type is used for the sandbox shortening benchmark. " *
                          "No slip boundary conditions are applied on bottom and left boundaries. Free " *
                          "slip boundary conditions are applied on top and right boundaries. Prescribed " *
                          "internal velocity boundary conditions are used to model a mobile wall. " *
                          "Boundary friction is also used along the left and bottom boundaries in " *
                          "addition to the interface between sand and mobile wall. Plastic failure is " *
                          "prevented from occurring in the mobile wall.",
                bools=Dict{Symbol, Bool}(
                    :use_boundary_friction_plasticity_model_sandbox => true, 
                    :no_yielding_in_mobile_wall => true
                    )
            ),
        option_ids[option_names.SandboxExtension] =>
            OptionState(
                option_name=string(option_names.SandboxExtension),
                description="This boundary condition type is used for the sandbox extension benchmark. " *
                          "No slip boundary conditions are applied on bottom and left boundaries. Free " *
                          "slip boundary conditions are applied on top and right boundaries. Prescribed " *
                          "internal velocity boundary conditions are used to model a mobile wall. " *
                          "Boundary friction is also used along the left and bottom boundaries in " *
                          "addition to the interface between sand and mobile wall. Plastic failure is " *
                          "prevented from occurring in the mobile wall and thin plate extending from " *
                          "mobile wall.",
                bools=Dict{Symbol, Bool}(
                    :use_boundary_friction_plasticity_model_sandbox => true, 
                    :no_yielding_in_mobile_wall => true,
                    :no_yielding_in_plate_extension => true
                    )
            ),
        option_ids[option_names.VerticalChannelFlowNonNewtonian] =>
            OptionState(
                option_name=string(option_names.VerticalChannelFlowNonNewtonian),
                description="This boundary condition type is used for vertical channel flow with " *
                          "prescribed pressure and velocity gradients along upper and lower boundaries. " *
                          "No-slip boundary condition are applied on side boundaries. Prescribed " *
                          "temperature is applied along top and bottom boundaries. Prescribed gradients " *
                          "are applied along side boundaries. Special regridding for is applied for " *
                          "channel flow. This bc type is used for non-newtonian channel flow benchmark " *
                          "with temperature along side boundaries set equal to zero.",
            ),
        option_ids[option_names.VerticalChannelFlowIsoviscous] =>
            OptionState(
                option_name=string(option_names.VerticalChannelFlowIsoviscous),
                description="This boundary conditions type is used for vertical channel flow with " *
                          "prescribed pressure and velocity gradients along upper and lower boundaries. " *
                          "No-slip boundary conditions are applied on side boundaries. Prescribed " *
                          "temperature along top and bottom boundaries. Prescribed gradients are applied " *
                          "along side boundaries. Special regridding for is applied for channel flow. " *
                          "This bc type is used for the isoviscous channel flow benchmark with thermal " *
                          "advection benchmark and temperature prescribed along side boundaries.",
            ),
        option_ids[option_names.VerticalChannelFlowShearHeating] =>
            OptionState(
                option_name=string(option_names.VerticalChannelFlowShearHeating),
                description="This boundary conditions type is used for vertical channel flow with " *
                          "prescribed pressure and velocity gradients along upper and lower boundaries. " *
                          "No-slip boundary conditions are applied on side boundaries. Prescribed " *
                          "temperature along top and bottom boundaries. Prescribed gradients are applied " *
                          "along side boundaries. Special regridding for is applied for channel flow. " *
                          "This boundary condition type is used for the channel flow benchmark with " *
                          "variable thermal conductivity and shear heating.",
            ),
        option_ids[option_names.FractureZoneContractionFixedBoundaries] =>
            OptionState(
                option_name=string(option_names.FractureZoneContractionFixedBoundaries),
                description="This boundary condition type is used for contraction models with fixed " *
                          "boundaries, symmetric contraction and fracture zone material setup. Free slip " *
                          "boundary conditions are applied along top, bottom and side boundaries. " *
                          "Orthogonal velocity is prescribed along side boundaries. Temperature is " *
                          "prescribed on top and bottom, and zero heat flux is applied along side " *
                          "boundaries.",
            ),
        option_ids[option_names.FractureZoneContractionFixedBoundariesAsymmetric] =>
            OptionState(
                option_name=string(option_names.FractureZoneContractionFixedBoundariesAsymmetric),
                description="This boundary condition type is used for contraction models with fixed " *
                          "boundaries, asymmetric contraction and fracture zone material setup. Free slip " *
                          "boundary conditions are applied along top, bottom and side boundaries. " *
                          "Orthogonal velocity is prescribed along side boundaries with zero inflow along " *
                          "the left boundary. Temperature is prescribed on top and bottom, and zero heat " *
                          "flux is applied along side boundaries.",
            ),
        option_ids[option_names.ViscoelasticContraction] =>
            OptionState(
                option_name=string(option_names.ViscoelasticContraction),
                description="This boundary condition type is used for contraction models with fixed " *
                          "boundaries, symmetric contraction and viscoelastic lithosphere material setup. " *
                          "Free slip boundary conditions are applied along top, bottom and side " *
                          "boundaries. Orthogonal velocity is prescribed along side boundaries. " *
                          "Temperature is prescribed on top and bottom and zero heat flux is applied " *
                          "along side boundaries.",
                bools=Dict{Symbol, Bool}(:stepping_velocity_bc => true)
            ),
        option_ids[option_names.ViscoelasticContractionAsymmetric] =>
            OptionState(
                option_name=string(option_names.ViscoelasticContractionAsymmetric),
                description="This boundary condition type is used for contraction models with fixed " *
                          "boundaries, asymmetric contraction and viscoelastic lithosphere material setup. " *
                          "Free slip boundary conditions are applied along top, bottom and side " *
                          "boundaries. Orthogonal velocity is prescribed along side boundaries with zero " *
                          "inflow along the left boundary. " *
                          "Temperature is prescribed on top and bottom and zero heat flux is applied " *
                          "along side boundaries.",
                bools=Dict{Symbol, Bool}(:stepping_velocity_bc => true)
            )
    )
end

end # module 