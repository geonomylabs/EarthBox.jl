module FlattenDensity

import EarthBox.ModelDataContainer: ModelData

"""
    flatten_staggered_grid_density_in_sticky_air_layer!(model::ModelData)::Nothing

Flatten grid density distribution for "air/water" boundary.

Density defined on basic grid is forced to conform to the density of
sticky air or sticky water. This avoids perturbations in the sticky
air/water weak layer.

Grid nodes are identified as sticky air or sticky water based on the
topography model being activated, the density being less than or equal to
1100 kg/m^3, radiogenic heat production being zero (except for vy grid),
the node being above or below the water surface and there being defined
sticky_air and sticky_water material types. If these conditions are not met
no flattening occurs.

# Arguments
- `model::ModelData`: The model data structure containing all necessary arrays and parameters

# Updated Arrays
## Updated arrays from group `model.stokes_continuity.arrays.density`:
- `rho1.array`: Density (kg/m^3) interpolated from markers to basic grid
- `rho1_vy.array`: Density (kg/m^3) interpolated from markers to vy grid
"""
function flatten_staggered_grid_density_in_sticky_air_layer!(model::ModelData)::Nothing
    # Get sticky air material ID
    sticky_air_ids = model.materials.dicts.matid_types["StickyAir"]
    matid_sticky_air = length(sticky_air_ids) > 0 ? first(sticky_air_ids) : -1

    # Get sticky water material ID
    sticky_water_ids = model.materials.dicts.matid_types["StickyWater"]
    matid_sticky_water = length(sticky_water_ids) > 0 ? first(sticky_water_ids) : -1

    if model.topography.parameters.model_options.iuse_topo.value == 1
        ynum = model.grids.parameters.geometry.ynum.value
        xnum = model.grids.parameters.geometry.xnum.value
        gridy_b = model.grids.arrays.basic.gridy_b.array
        hr1 = model.heat_equation.arrays.radiogenic_production.hr1.array
        mat_rho = model.materials.arrays.mat_rho.array
        y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
        
        # Density on basic grid
        rho1 = model.stokes_continuity.arrays.density.rho1.array
        for j in 1:xnum
            for i in 1:ynum
                if rho1[i,j] ≤ 1100.0 && gridy_b[i] < y_sealevel && hr1[i,j] == 0.0
                    # Only update density if sticky air has been defined
                    if matid_sticky_air > -1
                        rho1[i,j] = mat_rho[matid_sticky_air, 1]
                    end
                end

                if rho1[i,j] ≤ 1100.0 && gridy_b[i] ≥ y_sealevel && hr1[i,j] == 0.0
                    # Only update density if sticky water has been defined
                    if matid_sticky_water > -1
                        rho1[i,j] = mat_rho[matid_sticky_water, 1]
                    end
                end
            end
        end

        # Density on vy grid
        # Radiogenic heat production is not used for Vy density flattening
        # since it is not defined on the Vy grid
        rho1_vy = model.stokes_continuity.arrays.density.rho1_vy.array
        for j in 1:xnum+1
            for i in 1:ynum
                if rho1_vy[i,j] ≤ 1100.0 && gridy_b[i] < y_sealevel
                    # Only update density if sticky air has been defined
                    if matid_sticky_air > -1
                        rho1_vy[i,j] = mat_rho[matid_sticky_air, 1]
                    end
                end

                if rho1_vy[i,j] ≤ 1100.0 && gridy_b[i] ≥ y_sealevel
                    # Only update density if sticky water has been defined
                    if matid_sticky_water > -1
                        rho1_vy[i,j] = mat_rho[matid_sticky_water, 1]
                    end
                end
            end
        end
    end
    return nothing
end

end # module 