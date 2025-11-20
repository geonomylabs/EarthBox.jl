module Grids3dContainer

include("parameters/ParameterCollection.jl")
include("arrays/ArrayCollection.jl")

import EarthBox.EarthBoxDtypes: CollectionContainer
import .ArrayCollection: Arrays
import .ParameterCollection: Parameters

"""
    Grids3d

Data structure containing all 3d grid-related components for the EarthBox model.
This is a pure data container that holds the state of the grid system.
All operations on this data should be performed through the StaggeredGrid module.

# Fields
- `parameters`: Grid configuration parameters including options, geometry, and
  refinement settings
- `arrays`: Grid arrays containing coordinates and spacing for basic and staggered
  grids

# Composition

## Parameter Objects
- parameter object attributes: value, name, units and description

### Grid Option Parameters
- `parameters.grid_options.itype_grid`

### Geometry Parameters
- `parameters.geometry.xnum`
- `parameters.geometry.ynum`
- `parameters.geometry.znum`
- `parameters.geometry.xsize`
- `parameters.geometry.ysize`
- `parameters.geometry.zsize`
- `parameters.geometry.xstpavg`
- `parameters.geometry.ystpavg`
- `parameters.geometry.zstpavg`
- `parameters.geometry.ymin`
- `parameters.geometry.ymax`
- `parameters.geometry.xmin`
- `parameters.geometry.xmax`
- `parameters.geometry.zmin`
- `parameters.geometry.zmax`
- `parameters.geometry.xsize_start`

### Refinement Parameters
- `parameters.refinement.dx_highres`
- `parameters.refinement.xo_highres`
- `parameters.refinement.ixo_highres`
- `parameters.refinement.xf_highres`
- `parameters.refinement.dz_highres`
- `parameters.refinement.zo_highres`
- `parameters.refinement.izo_highres`
- `parameters.refinement.zf_highres`
- `parameters.refinement.dy_highres`
- `parameters.refinement.yf_highres`
- `parameters.refinement.iuse_trench`
- `parameters.refinement.trench_location`
- `parameters.refinement.iuse_refinement_delay`
- `parameters.refinement.refinement_time`
- `parameters.refinement.refinement_flag`
- `parameters.refinement.iuse_refinement_gap`
- `parameters.refinement.refinement_gap_start_time`
- `parameters.refinement.refinement_gap_end_time`

## Array Objects
- array object attributes: array, name, units and description)

### Basic Grid Arrays
- `arrays.basic.gridx_b`
- `arrays.basic.xstp_b`
- `arrays.basic.gridy_b`
- `arrays.basic.ystp_b`
- `arrays.basic.gridz_b`
- `arrays.basic.zstp_b`

### Staggered Vy Grid Arrays
- `arrays.staggered_vy.gridy_vy`
- `arrays.staggered_vy.ystp_vy`
- `arrays.staggered_vy.gridx_vy`
- `arrays.staggered_vy.xstp_vy`
- `arrays.staggered_vy.gridz_vy`
- `arrays.staggered_vy.zstp_vy`

### Staggered Vx Grid Arrays
- `arrays.staggered_vx.gridx_vx`
- `arrays.staggered_vx.xstp_vx`
- `arrays.staggered_vx.gridy_vx`
- `arrays.staggered_vx.ystp_vx`
- `arrays.staggered_vx.gridz_vx`
- `arrays.staggered_vx.zstp_vx`

### Staggered Vz Grid Arrays
- `arrays.staggered_vz.gridz_vz`
- `arrays.staggered_vz.zstp_vz`
- `arrays.staggered_vz.gridx_vz`
- `arrays.staggered_vz.xstp_vz`
- `arrays.staggered_vz.gridy_vz`
- `arrays.staggered_vz.ystp_vz`

### Pressure Grid Arrays
- `arrays.pressure.gridy_pr`
- `arrays.pressure.ystp_pr`
- `arrays.pressure.gridx_pr`
- `arrays.pressure.xstp_pr`
- `arrays.pressure.gridz_pr`
- `arrays.pressure.zstp_pr`

### Staggered Sxy Grid Arrays
- `arrays.staggered_sxy.gridx_sxy`
- `arrays.staggered_sxy.xstp_sxy`
- `arrays.staggered_sxy.gridy_sxy`
- `arrays.staggered_sxy.ystp_sxy`
- `arrays.staggered_sxy.gridz_sxy`
- `arrays.staggered_sxy.zstp_sxy`

### Staggered Sxz Grid Arrays
- `arrays.staggered_sxz.gridx_sxz`
- `arrays.staggered_sxz.xstp_sxz`
- `arrays.staggered_sxz.gridy_sxz`
- `arrays.staggered_sxz.ystp_sxz`
- `arrays.staggered_sxz.gridz_sxz`
- `arrays.staggered_sxz.zstp_sxz`

### Staggered Syz Grid Arrays
- `arrays.staggered_syz.gridx_syz`
- `arrays.staggered_syz.xstp_syz`
- `arrays.staggered_syz.gridy_syz`
- `arrays.staggered_syz.ystp_syz`
- `arrays.staggered_syz.gridz_syz`
- `arrays.staggered_syz.zstp_syz`

"""
mutable struct Grids3d <: CollectionContainer
    parameters::Parameters
    arrays::Arrays
end

function Grids3d(;
    ynum::Int64, 
    xnum::Int64, 
    znum::Int64, 
    ysize::Float64, 
    xsize::Float64, 
    zsize::Float64
)::Grids3d
    @assert ynum > 2 "Number of grid points in y-direction must be greater than 2"
    @assert xnum > 2 "Number of grid points in x-direction must be greater than 2"
    @assert znum > 2 "Number of grid points in z-direction must be greater than 2"
    @assert ysize > 0.0 "Domain size in y-direction must be positive"
    @assert xsize > 0.0 "Domain size in x-direction must be positive"
    @assert zsize > 0.0 "Domain size in z-direction must be positive"
    @assert ysize/ynum > 0.0 "Grid spacing in y-direction must be positive"
    @assert xsize/xnum > 0.0 "Grid spacing in x-direction must be positive"
    @assert zsize/znum > 0.0 "Grid spacing in z-direction must be positive"
    return Grids3d(
      Parameters(ynum, xnum, znum, ysize, xsize, zsize),
      Arrays(ynum, xnum, znum)
    )
end

end # module 