module ArrayLookupManager

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Arrays.ArrayTypes: ScalarArray2D
import ..DataNames: ScalarH5Datanames

struct ArrayLookup
    datanames::ScalarH5Datanames
    lookup::Dict{String, ScalarArray2D}
end

function ArrayLookup(model::ModelData)
    datanames = ScalarH5Datanames()
    lookup = build_lookup(model)
    return ArrayLookup(datanames, lookup)
end

function build_lookup(model::ModelData)::Dict{String, ScalarArray2D}
    lookup = Dict{String, ScalarArray2D}(
        datanames.TempC => model.heat_equation.arrays.temperature.tk1,
        datanames.log_eta_Pas => model.stokes_continuity.arrays.viscosity.etas1,
        datanames.Eii_log => model.stokes_continuity.arrays.strain_rate_and_spin.eii,
        datanames.pressure_GPa => model.stokes_continuity.arrays.pressure.pr1,
        datanames.Sxx_MPa => model.stokes_continuity.arrays.stress.sxx2,
        datanames.Sxy_MPa => model.stokes_continuity.arrays.stress.sxy2,
        datanames.plastics => model.stokes_continuity.arrays.plastic_def.plastics,
        datanames.plasticn => model.stokes_continuity.arrays.plastic_def.plasticn,
        datanames.velx_cmyr => model.stokes_continuity.arrays.basic_grid_velocity.vxb,
        datanames.vely_cmyr => model.stokes_continuity.arrays.basic_grid_velocity.vyb
    )
    return lookup
end

function get_array_info(
    lookup::ArrayLookup,
    dataname::String
)::Tuple{Matrix{Float64}, String, String}
    model_array = lookup.lookup[dataname]
    ebarray2d = model_array.array
    units = model_array.units
    grid_type = model_array.grid_type
    return ebarray2d, units, grid_type
end

end # module ArrayLookup 