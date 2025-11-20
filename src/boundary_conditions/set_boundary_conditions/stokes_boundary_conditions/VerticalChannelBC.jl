"""
    Boundary condition coefficients for infinite open vertical channel flow.
"""
module VerticalChannelBC

using EarthBox.ModelDataContainer: ModelData
using ..StokesGradientBC
using ..PressureBC

"""
    set_open_vertical_channel!(model::ModelData)::Nothing

Set boundary condition coefficient arrays for vertical channel flow.

The channel for this case is infinite and open at the top and bottom boundaries.
A vertical pressure gradient can be prescribed using the input pressure_bc
parameter, which is applied at the top boundary while zero pressure is applied
at the bottom boundary.
"""
function set_open_vertical_channel!(model::ModelData)::Nothing
    # Set pressure mode to prescribed pressure along top and bottom boundaries.
    PressureBC.set_pressure_bc_mode!(model, 1)
    # Set velocity boundary conditions along top boundary
    StokesGradientBC.set_zero_vertical_vx_gradient_at_top_boundary!(model)
    StokesGradientBC.set_zero_vertical_vy_gradient_in_top_cells!(model)
    # Set velocity boundary conditions along bottom boundary
    StokesGradientBC.set_zero_vertical_vx_gradient_at_bottom_boundary!(model)
    StokesGradientBC.set_zero_vertical_vy_gradient_at_bottom_cells!(model)
    return nothing
end

end # module
