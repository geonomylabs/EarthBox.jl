module StickyThickness

import EarthBox.ModelDataContainer: ModelData

function get_sticky_thickness(
    model::ModelData
)::Tuple{Float64, Float64}
    iuse_topo = model.topography.parameters.model_options.iuse_topo.value
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value

    if iuse_topo == 0 || ntimestep == 0
        sticky_air_thickness = model.geometry.parameters.sticky_air_geometry.thick_air.value
        sticky_thickness_left = sticky_air_thickness
        sticky_thickness_right = sticky_air_thickness
    else
        gridt = model.topography.arrays.gridt.array
        sticky_thickness_left = gridt[2, 1]  # Julia uses 1-based indexing
        sticky_thickness_right = gridt[2, end]
    end

    return (sticky_thickness_left, sticky_thickness_right)
end

end # module StickyThickness 