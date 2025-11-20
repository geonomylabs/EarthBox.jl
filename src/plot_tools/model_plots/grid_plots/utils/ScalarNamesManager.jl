module ScalarNamesManager

Base.@kwdef struct ScalarNames
    temperature::String = "temperature"
    viscosity::String = "viscosity"
    strainrate::String = "strainrate"
    pressure::String = "pressure"
    normal_stress::String = "normal_stress"
    shear_stress::String = "shear_stress"
    shear_plastic_failure::String = "shear_plastic_failure"
    normal_plastic_failure::String = "normal_plastic_failure"
    velocity_x::String = "velocity_x"
    velocity_y::String = "velocity_y"
    velocity_mag::String = "velocity_mag"
    density::String = "density"
    thermal_conductivity::String = "thermal_conductivity"
end

function get_list(sn::ScalarNames)::Vector{String}
    return [getfield(sn, field) for field in fieldnames(ScalarNames)]
end

function get_default_scalar_plot_parameters()
    sn = ScalarNames()
    return Dict(
        # Assuming Celsius
        sn.temperature => (
            contour_interval = 50.0,
            minimum_value = 0.01,
            maximum_value = 1300.0,
            plot_contours = true,
            grid_plot_type = "nomesh"
        ),
        # Assuming log10(Pa.s)
        sn.viscosity => (
            contour_interval = 0.5,
            minimum_value = 18.0,
            maximum_value = 26.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming log10(1/s)
        sn.strainrate => (
            contour_interval = 0.5,
            minimum_value = -18.0,
            maximum_value = -14.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming GPa
        sn.pressure => (
            contour_interval = 1.0,
            minimum_value = 0.0,
            maximum_value = 10.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming MPa
        sn.normal_stress => (
            contour_interval = 5.0,
            minimum_value = -10.0,
            maximum_value = 10.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming MPa
        sn.shear_stress => (
            contour_interval = 5.0,
            minimum_value = -10.0,
            maximum_value = 10.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming 0/1
        sn.shear_plastic_failure => (
            contour_interval = 0.1,
            minimum_value = 0.0,
            maximum_value = 1.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming 0/1
        sn.normal_plastic_failure => (
            contour_interval = 0.1,
            minimum_value = 0.0,
            maximum_value = 1.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming cm/yr
        sn.velocity_x => (
            contour_interval = 0.25,
            minimum_value = -2.5,
            maximum_value = 2.5,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming cm/yr
        sn.velocity_y => (
            contour_interval = 0.25,
            minimum_value = -2.5,
            maximum_value = 2.5,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming cm/yr
        sn.velocity_mag => (
            contour_interval = 0.25,
            minimum_value = 0.0,
            maximum_value = 5.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming kg/m^3
        sn.density => (
            contour_interval = 50.0,
            minimum_value = 2400.0,
            maximum_value = 3400.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        ),
        # Assuming W/m/K
        sn.thermal_conductivity => (
            contour_interval = 0.5,
            minimum_value = 0.0,
            maximum_value = 20.0,
            plot_contours = false,
            grid_plot_type = "nomesh"
        )
    )
end

end # module 