# EarthBox Installation

## Operating System

Although the core functionality of Julia and EarthBox will work on Windows
and MacOS operating systems, we recommend using a Linux operating system since 
Linux yields better performance for scientific software, and the installation 
of complex scientific software often involves less overhead.

EarthBox was developed and tested on Debian-based Linux distributions like 
Ubuntu and PopOS. Therefore, you will see a strong bias toward Debian-based package 
management and the Bash shell in the EarthBox documentation. We recommend using Ubuntu
or PopOS since EarthBox was extensively tested on these operating systems. 

If you are working on Windows we recommend using VirtualBox to build a virtual Linux
machine. See [VirtualBox](https://www.virtualbox.org).


## Julia Installation

The easiest way to get started with Julia is to use `juliaup` which involves
executing the following command on Linux systems:

```bash
curl -fsSL https://install.julialang.org | sh
```

Juliaup will automatically customize your `~/.bashrc` file for using Julia.

See [Installing Julia](https://julialang.org/install/) for more details. 

Alternatively, you can install Julia directly from the latest tarball. For example:

```bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.1-linux-x86_64.tar.gz
tar zxvf julia-1.12.1-linux-x86_64.tar.gz
```

See [Julia Tarball](https://julialang.org/downloads/platform/#linux_and_freebsd) for
more details. Custom installations of Julia will require adding the path to the Julia 
bin in your `~/.bashrc` file:

``` bash
export PATH="path/to/the/unpacked/julia/bin:$PATH"
```

and then sourcing your `~/.bashrc`:

```bash
source ~/.bashrc
```

It is often useful to define a symbolic link to the Julia directory, which includes 
the version number in the directory name, to eliminate the need to update your 
`~/.bashrc` with when a new version is installed. For example,

```bash
ln -sfn julia-1.11.2 julia
```

You can customize your Julia setup by defining the following environment variables
in your `~/.bashrc`:

```bash
# Tell Julia to use all available threads
export JULIA_NUM_THREADS=auto
# Define where packages, artifacts etc... will be stored
export JULIA_DEPOT_PATH="/path/to/your/custom/julia-depot"
export JULIA_SCRATCH="/path/to/your/custom/julia-scratch"                   
export JULIA_ARTIFACTS_DIR="/path/to/your/custom/julia-artifacts"
```

When using MUMPS it is often a good idea to restrict the number of OMP threads
to 1 in your `~/.bashrc` especially if you are using a locally compiled version:

```bash
export OMP_NUM_THREADS=1
```

## Basic EarthBox Installation

You can install the registered version of EarthBox via the Julia REPL:

```julia
julia> ]
(@v1.x) pkg> add EarthBox
```


## Installing MPI in Base Julia Environment

If you plan on using the MUMPS solver you will need to install MPI in your base 
Julia environment, configure MPI using MPIPreferences and install `mpiexecjl`:

[1] Add MPI to your base environment:

```julia
julia> ]
(@v1.12) pkg> add MPI@0.20.22
```

[2] And MPIPreferences to your base environment:

```julia
julia> ]
(@v1.12) pkg> add MPIPreferences@0.1.11
```

!!! warning "version matching with EarthBox"
    It is a good idea to match the version of MPI and MPIPreferences used by
    EarthBox.

See [Stable MPI Configuration](https://juliaparallel.org/MPI.jl/stable/usage/)
for more details.

[3] configure MPI to use the Julia OpenMPI binary:

```julia
julia> using MPIPreferences
julia> MPIPreferences.use_jll_binary("OpenMPI_jll")
```

!!! tip "error message"
    An error message is generated instructing you to restart Julia. This
    is OK, the task was completed properly. Just restart the Julia REPL.

[4] Install `mpiexecjl` using the following:

```julia
julia> using MPI
julia> MPI.install_mpiexecjl()
```

!!! tip "note the path to mpiexecjl"
    make a note of the path to the directory containing `mpiexecjl` produce as screen 
    output since you will add this path to your `.bashrc`

[5] Add the path to `mpiexecjl` to your `~/.bashrc`

```bash
export PATH="/path/to/mpiexecjl/bin:$PATH"
```

```bash
source ~/.bashrc
```


## EarthBox Installation for Development and Experimentation

The  following instructions are for installing the EarthBox project for local 
experimentation that involves modifying the code.

[1] Download the EarthBox project directory. 

[2] Navigate to the EarthBox directory and ensure that `Manifest.toml` and 
`LocalPreferences.toml` files are deleted. Then enter the REPL with EarthBox as 
the active project:

```bash
julia --project=.
```

[3] From inside the REPL instantiate the EarthBox project:

```
julia> using Pkg
julia> Pkg.instantiate()
```

This may take some time to complete especially if you are starting with a fresh
depot. 

[4] With the active EarthBox project, configure MPI to use Julia binaries for OpenMP.

```julia
julia> using MPIPreferences
julia> MPIPreferences.use_jll_binary("OpenMPI_jll")
```

MPI.jl should also be installed in the base Julia environment (see 
[Installing MPI in Base Julia Environment](@ref)). These MPI configuration steps
are required for the use of MUMPS.jl with EarthBox.

[5] Add EarthBox to your current default Julia environment for development.
First exit out of the REPL if you are in it (you can use `exit()` from the REPL). 
Then re-enter the REPL with the base Julia environment and add the EarthBox 
project as a development project:

```julia
julia> ]
(@v1.x) pkg> dev /path/to/EarthBox
```

This will allow you to access EarthBox via

``` julia
using EarthBox
```

After these steps are completed you will have a basic working version of EarthBox
that used Julia binaries for MPI and MUMPS.

!!! warning "Parallel Analysis with MUMPS Using PT-Scotch and ParMetis"
    Parallel analysis with PT-Scotch and ParMetis for the MUMPS solver may not 
    work on your system when using JUlia binaries. You may need to set the MUMPS 
    `analysis_method` to SERIAL to work around these issues. However, having parallel 
    analysis can significantly increase MUMPS performance especially with large 
    systems of equations so it is worth the effort to get these features working. 
    See the [Configuring MPI.jl to Use System Binaries](@ref) and 
    [Installing Custom MUMPS Libraries](@ref) for guidance on using system MPI 
    and a local compiled version of MUMPS on Debian-type systems (e.g. Ubuntu, PopOS).