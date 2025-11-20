module DataStructures

""" Sediment transport parameters.

# Fields
- `subaerial_slope_diffusivity::Float64`: subaerial slope diffusivity (m^2/s)
- `precipitation_rate::Float64`: precipitation rate (m/s)
- `subaerial_transport_coefficient::Float64`: subaerial transport coefficient
- `submarine_slope_diffusivity::Float64`: submarine slope diffusivity (m^2/s)
- `submarine_diffusion_decay_depth::Float64`: submarine diffusion decay depth (m)
- `transport_timestep::Float64`: transport timestep (s)
- `transport_model_duration::Float64`: transport model duration (s)
- `porosity_initial::Float64`: initial porosity (fraction)
- `depth_decay_term::Float64`: depth decay term (1/m)
"""
struct SedimentTransportParameters
    subaerial_slope_diffusivity::Float64  # m^2/s
    precipitation_rate::Float64           # m/s
    subaerial_transport_coefficient::Float64
    submarine_slope_diffusivity::Float64  # m^2/s
    submarine_diffusion_decay_depth::Float64  # m
    transport_timestep::Float64           # s
    transport_model_duration::Float64     # s
    porosity_initial::Float64             # fraction
    depth_decay_term::Float64             # 1/m
end

end # module 