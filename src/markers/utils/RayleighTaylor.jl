module RayleighTaylor

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData

using Printf

function perturb_marker_y!(model::ModelData)
    xsize = model.grids.parameters.geometry.xsize.value
    wavelength = model.geometry.parameters.rayleigh_taylor.wave_length_lambda.value
    amplitude_initial = model.geometry.parameters.rayleigh_taylor.amplitude_initial.value

    print_rayleigh_taylor_info(xsize, wavelength, amplitude_initial)

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marknum = model.markers.parameters.distribution.marknum.value

    for imarker in 1:marknum
        mx_coor = marker_x[imarker]
        my_coor = marker_y[imarker]
        marker_y[imarker] = rayleigh_taylor_instability_ycoordinate(
            mx_coor, my_coor, xsize, wavelength, amplitude_initial
        )
    end
end

function print_rayleigh_taylor_info(xsize, wavelength, amplitude_initial)
    print_info("", level=1)
    print_info("Adjusting coordinates for Rayleigh-Taylor Instability", level=1)
    print_info("xsize: $xsize", level=2)
    print_info("wavelength: $wavelength", level=2)
    print_info("amplitude_initial: $amplitude_initial", level=2)
    print_info("", level=1)
end

""" Update marker y-coordinate for Rayleigh-Taylor instability benchmark.

# Arguments
- `mx_coor`: X-coordinate of marker
- `my_coor`: Y-coordinate of marker
- `xsize`: Size of domain in x-direction
- `wavelength`: Wavelength of perturbation
- `amplitude_initial`: Initial amplitude of perturbation

# Returns
- `my_coor_new`: Updated y-coordinate
"""
function rayleigh_taylor_instability_ycoordinate(
    mx_coor::Float64,
    my_coor::Float64,
    xsize::Float64,
    wavelength::Float64,
    amplitude_initial::Float64
)::Float64
    my_coor_new = (
        my_coor 
        - cos(2.0 * π * (mx_coor - xsize/2.0) / wavelength) * amplitude_initial
        )
    return my_coor_new
end

end # module RayleighTaylor 