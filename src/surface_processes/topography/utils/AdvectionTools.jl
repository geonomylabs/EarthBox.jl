module AdvectionTools

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Interpolation: GridToMarker
import EarthBox.Interpolation.MarkerGridMapping: get_index_of_left_node
import EarthBox.Interpolation.MarkerGridMapping: upr_left_x_mapping_vx_grid
import EarthBox.Interpolation.MarkerGridMapping: upr_left_y_mapping_vx_grid
import EarthBox.Interpolation.MarkerGridMapping: upr_left_x_mapping_vy_grid
import EarthBox.Interpolation.MarkerGridMapping: upr_left_y_mapping_vy_grid
import EarthBox: ConversionFuncs
import EarthBox.SedimentThickness: calculate_sediment_thickness_from_markers
import EarthBox.Compaction.ApplyCompaction: apply_compaction_model
import EarthBox.DataStructures: SedimentTransportParameters

function set_topo_velocity_to_zeros!(toponum::Int, gridt::Array{Float64,2})::Nothing
    for i in 1:toponum
        gridt[4, i] = 0.0
        gridt[5, i] = 0.0
        gridt[6, i] = 0.0
    end
    return nothing
end

""" Calculate advection velocity at topography points.

# Arguments
- `model::ModelData`: Model data container
- `gridt::Array{Float64,2}`: Multi-dimensional topography grid array with elevation, velocity and
  antidiffusion information

# Updated Arrays
- `gridt[4, toponum]`: x-component of velocity
- `gridt[5, toponum]`: y-component of velocity

# Interpolation Method
```
  xn    V(xn,yn)--------------------V(xn+1,yn)
          ?           ^                  ?
          ?           ?                  ?
          ?          dy                  ?
          ?           ?                  ?
          ?           v                  ?
          ?<----dx--->o Mrho(xm,ym)       ?
          ?                              ?
          ?                              ?
  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)
```
"""
function calc_advection_veloc_from_basic_grid!(
    model::ModelData,
    gridt::Matrix{Float64}
)::Nothing
    dy_upper = 0.0

    toponum = model.topography.parameters.topo_grid.toponum.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array

    for i in 1:toponum
        gridt[4, i] = 0.0  # x-component of velocity
        gridt[5, i] = 0.0  # y-component of velocity
        gridt[6, i] = 0.0  # antidiffusion correction
    end

    for i in 1:toponum
        if topo_node_inside_grid(model, i, gridt)
            xcur = gridt[1, i]
            ycur = gridt[2, i] + dy_upper
            x_index_upr_left_basic = get_index_of_left_node(xcur, gridx_b, xnum)
            y_index_upr_left_basic = get_index_of_left_node(ycur, gridy_b, ynum)
            # X-component of velocity
            gridt[4, i] = get_velocity_x(
                model,
                xcur, ycur,
                x_index_upr_left_basic,
                y_index_upr_left_basic
            )
            # Y-component of velocity
            gridt[5, i] = get_velocity_y(
                model,
                xcur, ycur,
                x_index_upr_left_basic,
                y_index_upr_left_basic
            )
        end
    end
    return nothing
end

function get_velocity_x(
    model::ModelData,
    xcur::Float64,
    ycur::Float64,
    x_index_upr_left_basic::Int32,
    y_index_upr_left_basic::Int32
)::Float64
    x_index_upr_left, dx_upr_left = upr_left_x_mapping_vx_grid(
        x_index_upr_left_basic,
        xcur,
        model.grids.arrays.basic.gridx_b.array,
        model.grids.arrays.basic.xstp_b.array
    )

    y_index_upr_left, dy_upr_left = upr_left_y_mapping_vx_grid(
        y_index_upr_left_basic,
        ycur,
        model.grids.parameters.geometry.ynum.value,
        model.grids.arrays.staggered_vx.gridy_vx.array,
        model.grids.arrays.staggered_vx.ystp_vx.array
    )

    velocity_x = GridToMarker.get_marker_value(
        y_index_upr_left,
        x_index_upr_left,
        dy_upr_left,
        dx_upr_left,
        model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    )

    return velocity_x
end

function get_velocity_y(
    model::ModelData,
    xcur::Float64,
    ycur::Float64,
    x_index_upr_left_basic::Int32,
    y_index_upr_left_basic::Int32
)::Float64
    x_index_upr_left, dx_upr_left = upr_left_x_mapping_vy_grid(
        x_index_upr_left_basic,
        xcur,
        model.grids.parameters.geometry.xnum.value,
        model.grids.arrays.staggered_vy.gridx_vy.array,
        model.grids.arrays.staggered_vy.xstp_vy.array
    )
    
    y_index_upr_left, dy_upr_left = upr_left_y_mapping_vy_grid(
        y_index_upr_left_basic,
        ycur,
        model.grids.arrays.basic.gridy_b.array,
        model.grids.arrays.basic.ystp_b.array
    )

    velocity_y = GridToMarker.get_marker_value(
        y_index_upr_left,
        x_index_upr_left,
        dy_upr_left,
        dx_upr_left,
        model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array
    )

    return velocity_y
end

function advect_vertically!(
    model::ModelData;
    use_stokes::Bool=false,
    use_uniform_erosion_and_deposition::Bool=false
)::Nothing
    if use_stokes
        advect_vertically_using_stokes_velocity!(model)
    end
    if use_uniform_erosion_and_deposition
        advect_vertically_using_uniform_erosion_and_deposition!(model)
    end
    return nothing
end

""" Advect topography vertically using Stokes y-velocity.

# Updated Arrays
- `gridt[2, toponum]`: Topography grid array with y-coordinate of topography nodes (elevation) (m)
"""
function advect_vertically_using_stokes_velocity!(model::ModelData)::Nothing
    gridt = get_topography_array(model)
    timestep = get_model_time_step(model)
    toponum = size(gridt, 2)
    
    for i in 1:toponum
        velocity_topo_y_stokes = gridt[5, i]
        y_initial_topo = gridt[2, i]
        gridt[2, i] = calculate_updated_topography_node_y_coordinate(
            y_initial_topo, velocity_topo_y_stokes, timestep)
    end
    return nothing
end

""" Advect topography vertically using erosion, deposition rates.

This function uses the y-coordinate of water depth to determine if the
topography node is below water level. If the node is above water level,
the node is eroded at a uniform rate. If the node is below water level,
the node is deposited at a uniform rate.

# Updated Arrays
- `gridt[2, toponum]`: Topography grid array with elevation (m)
"""
function advect_vertically_using_uniform_erosion_and_deposition!(
    model::ModelData
)::Nothing
    erosion_rate, sedimentation_rate = get_erosion_and_sedimentation_rates(model)
    gridt = get_topography_array(model)
    y_sealevel = get_water_level_y_coordinate(model)
    timestep = get_model_time_step(model)
    toponum = size(gridt, 2)
    
    for i in 1:toponum
        y_topo_initial = gridt[2, i]
        velocity_y_topo = calculate_uniform_sedimentation_or_erosion_velocity(
            y_topo_initial, y_sealevel, erosion_rate, sedimentation_rate)
        gridt[2, i] = calculate_updated_topography_node_y_coordinate(
            y_topo_initial, velocity_y_topo, timestep)
    end
    return nothing
end

""" Advect topography vertically using extrusion velocity.

# Updated Arrays
- `gridt[2, toponum]`: Topography grid array with y-coordinate of topography nodes (elevation) (m)
"""
function advect_vertically_using_extrusion!(model::ModelData)::Nothing
    gridt = get_topography_array(model)
    timestep = get_model_time_step(model)
    toponum = size(gridt, 2)
    
    for i in 1:toponum
        y_topo_initial = gridt[2, i]
        extrusion_thickness = gridt[7, i]
        velocity_y_topo_extrusion = update_velocity_for_extrusion(
            extrusion_thickness, timestep, 0.0)
        gridt[2, i] = calculate_updated_topography_node_y_coordinate(
            y_topo_initial, velocity_y_topo_extrusion, timestep)
    end
    return nothing
end

""" Advect topography vertically using salt_deposition rate.

# Updated Arrays
- `gridt[2, toponum]`: Topography grid array with y-coordinate of topography nodes (elevation) (m)
"""
function advect_topography_vertically_for_salt_deposition!(model::ModelData)::Nothing
    salt_start_time, salt_end_time, salt_deposition_rate_m_s = get_salt_deposition_info(model)
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    
    println(">> Working on salt deposition at timesum (Myr): ",
            ConversionFuncs.seconds_to_years(timesum)/1e6)
    println(">> Start time (Myr): ",
            ConversionFuncs.seconds_to_years(salt_start_time)/1e6)
    println(">> End time (Myr): ",
            ConversionFuncs.seconds_to_years(salt_end_time)/1e6)
    println(">> Salt deposition rate (m/s): ", salt_deposition_rate_m_s)
    println(">> Water level (m): ", y_sealevel)

    if salt_start_time <= timesum <= salt_end_time
        gridt = get_topography_array(model)
        topo_gridx = copy(gridt[1, :])
        topo_gridy = copy(gridt[2, :])
        topo_gridy_initial = copy(gridt[2, :])

        timestep = get_model_time_step(model)
        toponum = size(gridt, 2)
        
        for i in 1:toponum
            y_topo_initial = gridt[2, i]
            topo_gridy_value = calculate_updated_topography_node_y_coordinate(
                y_topo_initial, -salt_deposition_rate_m_s, timestep)
            topo_gridy[i] = max(topo_gridy_value, y_sealevel)
        end

        sediment_thickness_initial = calculate_sediment_thickness_from_markers(
            model, topo_gridx)

        salt_decompaction_parameters = get_salt_decompaction_parameters()

        _, _, _ = apply_compaction_model(
            model, topo_gridx, topo_gridy, topo_gridy_initial,
            sediment_thickness_initial,
            salt_decompaction_parameters
        )

        for i in 1:toponum
            gridt[2, i] = topo_gridy[i]
        end
    end
    return nothing
end

""" Get sediment transport parameters.

Only porosity_initial and depth_decay_term are used in the
SedimentTransportParameters data structure for the salt model.

These parameters are set to have no decompaction in salt
"""
function get_salt_decompaction_parameters()::SedimentTransportParameters
    porosity_initial = 0.0  # fraction
    decay_depth_term = 1.0/2000.0  # 1/m
    return SedimentTransportParameters(
        subaerial_slope_diffusivity=0.0,
        precipitation_rate=0.0,
        subaerial_transport_coefficient=0.0,
        submarine_slope_diffusivity=0.0,
        submarine_diffusion_decay_depth=0.0,
        transport_timestep=0.0,
        transport_model_duration=0.0,
        porosity_initial=porosity_initial,
        depth_decay_term=decay_depth_term
    )
end

""" Get salt deposition information.

# Returns
- `salt_start_time::Float64`: Time when salt deposition starts (s)
- `salt_end_time::Float64`: Time when salt deposition ends (s)
- `salt_deposition_rate_m_s::Float64`: Salt deposition rate (m/s)
"""
function get_salt_deposition_info(
    model::ModelData
)::Tuple{Float64,Float64,Float64}
    salt_start_time = model.topography.parameters.depo_and_erosion_rates.salt_start_time.value
    salt_start_time = ConversionFuncs.millions_of_years_to_seconds(salt_start_time)
    salt_end_time = model.topography.parameters.depo_and_erosion_rates.salt_end_time.value
    salt_end_time = ConversionFuncs.millions_of_years_to_seconds(salt_end_time)
    salt_deposition_rate_m_s = model.topography.parameters.depo_and_erosion_rates.salt_deposition_rate.value
    return salt_start_time, salt_end_time, salt_deposition_rate_m_s
end

function get_model_time_step(model::ModelData)::Float64
    return model.timestep.parameters.main_time_loop.timestep.value
end

function get_water_level_y_coordinate(model::ModelData)::Float64
    return model.topography.parameters.sealevel.y_sealevel.value
end

function get_topography_array(model::ModelData)::Array{Float64,2}
    return model.topography.arrays.gridt.array
end

function get_erosion_and_sedimentation_rates(
    model::ModelData
)::Tuple{Float64,Float64}
    erosion_rate = model.topography.parameters.depo_and_erosion_rates.erosion_rate.value
    sedimentation_rate = model.topography.parameters.depo_and_erosion_rates.sedimentation_rate.value
    return erosion_rate, sedimentation_rate
end

""" Update velocity for extrusion.

# Arguments
- `thickness_of_extruded_material::Float64`: Thickness of extruded material (m)
- `timestep::Float64`: Timestep (s)
- `velocity_y_topo::Float64`: Velocity of topography node (m/s)

# Returns
- `velocity_y_topo::Float64`: Updated velocity of topography node (m/s)
"""
function update_velocity_for_extrusion(
    thickness_of_extruded_material::Float64,
    timestep::Float64,
    velocity_y_topo::Float64
)::Float64
    velocity_y_extrusion = calculate_extrusion_velocity(
        thickness_of_extruded_material, timestep)
    velocity_y_topo = velocity_y_topo - velocity_y_extrusion
    return velocity_y_topo
end

""" Calculate extrusion velocity.

# Returns
- `velocity_y_extrusion::Float64`: Velocity of extruded material (m/s)
"""
function calculate_extrusion_velocity(
    thickness_of_extruded_material::Float64,
    timestep::Float64
)::Float64
    velocity_y_extrusion = thickness_of_extruded_material/timestep
    return velocity_y_extrusion
end

""" Calculate updated topography node y-coordinate.

# Arguments
- `y_topo_initial::Float64`: Initial y-coordinate of topography node (m)
- `velocity_y_topo::Float64`: Velocity of topography node (m/s)
- `timestep::Float64`: Timestep (s)

# Returns
- `y_topo_updated::Float64`: Updated y-coordinate of topography node (m)
"""
function calculate_updated_topography_node_y_coordinate(
    y_topo_initial::Float64,
    velocity_y_topo::Float64,
    timestep::Float64
)::Float64
    y_topo_updated = y_topo_initial + velocity_y_topo*timestep
    return y_topo_updated
end

""" Update velocity using Stokes y-velocity.

# Arguments
- `velocity_y_topo::Float64`: Velocity of topography node (m/s)
- `velocity_y_stokes::Float64`: Y-component of Stokes velocity (m/s)

# Returns
- `velocity_y_topo::Float64`: Updated velocity of topography node (m/s)
"""
function update_velocity_using_stokes(
    velocity_y_topo::Float64,
    velocity_y_stokes::Float64
)::Float64
    velocity_y_topo = velocity_y_stokes + velocity_y_topo
    return velocity_y_topo
end

""" Calculate uniform sedimentation or erosion velocity.

If the topography node is below water level, the node is deposited at a
uniform rate. If the node is above water level, the node is eroded at a
uniform rate.

# Arguments
- `y_topography::Float64`: Y-coordinate of topography node (m)
- `y_water_level::Float64`: Y-coordinate of water level (m)
- `erosion_rate::Float64`: Erosion rate (m/s)
- `sedimentation_rate::Float64`: Sedimentation rate (m/s)

# Returns
- `velocity_y_topo::Float64`: Velocity of topography node (m/s)
"""
function calculate_uniform_sedimentation_or_erosion_velocity(
    y_topography::Float64,
    y_water_level::Float64,
    erosion_rate::Float64,
    sedimentation_rate::Float64
)::Float64
    if y_topography < y_water_level
        velocity_y_topo = erosion_rate
    else
        velocity_y_topo = -abs(sedimentation_rate)
    end
    return velocity_y_topo
end

""" Advect topography horizontally using anti-diffusion approach.

# Updated Arrays
- `gridt[2,toponum]`: Topography grid array with y-coordinate of topography nodes (elevation) (m)
- `gridt[3,toponum]`: Topography grid array with updated FCT elevation (m)
- `gridt[6,toponum]`: Antidiffusion correction (m)
"""
function advect_horizontally_with_antidiffusion!(model::ModelData)::Nothing
    dx_topo = model.topography.parameters.topo_grid.dx_topo.value
    toponum = model.topography.parameters.topo_grid.toponum.value
    gridt = model.topography.arrays.gridt.array
    dt, nt = define_topo_advection_timestep(model, gridt)
    
    # Defining FCT parameter MU
    mu = 1/8
    
    # Advect topography with FCT
    for _t in 1:nt
        # Step 0: Set new profile
        for i in 1:toponum
            gridt[3, i] = gridt[2, i]
        end
        
        # Step 1: Transport + numerical diffusion stage
        for i in 2:(toponum-1)
            # Defining FCT parameters EPS and NU
            eps = gridt[4, i]*dt/dx_topo  # Using Vx horizontal velocity
            nu = 1.0/8.0 + (eps^2.0)/2.0
            # Change topography
            gridt[3, i] = (
                gridt[2, i]
                - eps/2*(gridt[2, i+1] - gridt[2, i-1])
                + nu*(gridt[2, i+1] - 2*gridt[2, i] + gridt[2, i-1])
            )
        end
        
        # Step 2: Antidiffusion stage
        # Antidiffusion flow for the first cell
        gridt[6, 1] = 0.0
        for i in 2:(toponum-2)
            # Corrected antidiffusion flow for current cell
            delt0 = gridt[3, i] - gridt[3, i-1]
            delt1 = gridt[3, i+1] - gridt[3, i]
            delt2 = gridt[3, i+2] - gridt[3, i+1]
            s = sign(delt1)
            gridt[6, i] = s*max(0.0, min(s*delt2, s*delt0), mu*abs(delt1))
            gridt[2, i] = gridt[3, i] - gridt[6, i] + gridt[6, i-1]
        end
    end
    return nothing
end

""" Define topography advection timestep.

# Returns
- `dt::Float64`: Topography advection timestep
- `nt::Int`: Number of advection timesteps
"""
function define_topo_advection_timestep(
    model::ModelData,
    gridt::Array{Float64,2}
)::Tuple{Float64,Int}
    dx_topo = model.topography.parameters.topo_grid.dx_topo.value
    toponum = model.topography.parameters.topo_grid.toponum.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    
    # Define maximal horizontal velocity at topography nodes
    vxmax = 1e-32
    for i in 1:toponum
        vx = abs(gridt[4, i])
        if vx > vxmax
            vxmax = vx
        end
    end
    
    dt = timestep
    nt = 1
    if vxmax > 0
        dt = dx_topo/vxmax
        if dt < timestep
            nt = floor(Int, timestep/dt - 0.5) + 1
            dt = timestep/nt
        else
            dt = timestep
        end
    end
    return dt, nt
end

""" Determine if topography node is inside model domain.

# Returns
- `check::Bool`: True if node is inside model
"""
function topo_node_inside_grid(
    model::ModelData,
    i::Int,
    gridt::Array{Float64,2}
)::Bool
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array
    
    check = false
    if (gridx_b[xnum] >= gridt[1, i] >= gridx_b[1] &&
        gridy_b[ynum] >= gridt[2, i] >= gridy_b[1])
        check = true
    end
    return check
end

end # module 