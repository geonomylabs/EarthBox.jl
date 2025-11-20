module ExtensionFreeSlipAlongBottom

import EarthBox.ModelDataContainer: ModelData

"""
    calc_velocity!(model::ModelData)

Calculate velocity at top and bottom boundaries of extension models.

Boundary velocity is calculated to conserve mass based on extension velocity
along side boundaries.

# Updated Parameters
- `model.bcs.parameters.velocity`
  - `vyu.value`: Conservative vertical velocity (m/s) along upper boundary of
    extension model.
  - `vyl.value`: Conservative vertical velocity (m/s) along lower boundary of
    extension model.
"""
function calc_velocity!(model::ModelData)::Nothing
    velocity_extension_full = model.bcs.parameters.velocity.full_velocity_extension.value
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    model.bcs.parameters.velocity.vyu.value = velocity_extension_full * ysize / xsize
    model.bcs.parameters.velocity.vyl.value = 0.0
    return nothing
end

end # module 