module ScalingFactors

function calculate_residual_scaling_factors(
    density_reference::Float64,
    viscosity_reference::Float64,
    gravitation_acceleration::Float64,
    xstpavg::Float64
)::NamedTuple{(:stokesscale, :continscale), Tuple{Float64, Float64}}
    # Defining scale for Stokes residuals from y-Stokes equation
    # dSIGMAij/dj-dP/di=-RHO*gi=0  => Stokes scale=abs(RHO*gi)
    stokesscale = density_reference * gravitation_acceleration
    # Defining scale for Continuity residuals from y-Stokes equation
    # dvx/dx+dvy/dy=0 can be transformed to 2ETA(dvx/dx+dvy/dy)/dx=0 
    # which is similar to dSIGMAij/dj and has scale given above
    # therefore continuity scale = scale=abs(RHO*gi/ETA*dx)
    continscale = density_reference * gravitation_acceleration / viscosity_reference * xstpavg
    return stokesscale, continscale
end

end