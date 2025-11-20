# Configuring MPI.jl to Use System Binaries

Performance can often be improved by configuring `MPI.jl` to use system MPI 
libraries. This section provides guidance on how to do this on Debian Linux 
operating systems (e.g. Ubuntu/PopOS) assuming that `MPI.jl` and `MPIPreferences.jl`
have already been installed for EarthBox and the base Julia environment (see 
[EarthBox Installation](@ref)): 

[1] Install required System libraries (OpenMPI)

A system-level OpenMPI MPI implementation can be installed on Debian Linux systems 
(Ubuntu/PopOS systems) using the following command:

```bash
sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev libopenmpi3
```

[2] Use `MPIPreferences` to configure `MPI.jl` to use the system binary:

If you are working with a development version of the EarthBox project, 
navigate to the project directory and run the REPL with the project activated:

```bash
julia --project=.
```

Now change the MPI configuration to use system binaries:

```julia
using MPIPreferences
MPIPreferences.use_system_binary()
```

Which will generate output that looks like the following:

```julia
┌ Info: MPI implementation identified
│   libmpi = "libmpi"
│   version_string = "Open MPI v4.1.2, package: Debian OpenMPI, ident: 4.1.2, repo rev: v4.1.2, Nov 24, 2021\0"
│   impl = "OpenMPI"
│   version = v"4.1.2"
└   abi = "OpenMPI"
┌ Info: MPIPreferences changed
│   binary = "system"
│   libmpi = "libmpi"
│   abi = "OpenMPI"
│   mpiexec = "mpiexec"
│   preloads = Any[]
└   preloads_env_switch = nothing
```

!!! tip "going back to JLL binaries"
    MPIPreferences.use_jll_binary("OpenMPI_jll")
 
[3] Re-instantiate your project:

```julia
Pkg.resolve()
Pkg.instantiate()
```
