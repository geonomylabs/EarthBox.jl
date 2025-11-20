module RunTools

import EarthBox: GetPaths
import EarthBox.GetArgs: get_case_name_from_cl_args
import EarthBox.GetArgs: get_root_path_from_args
import EarthBox.GetArgs: get_runit_actions_from_cl_args
import EarthBox.GetArgs: get_model_output_path_from_cl_args
import EarthBox.GetArgs: manage_earthbox_project_path
import EarthBox.GetArgs: get_istart, get_iend

""" 
    run_earthbox(;
        model_dir::String,
        eb_project_path_from_script::Union{String, Nothing} = nothing,
        root_path_output_from_script::Union{String, Nothing} = nothing
    )

Run the model script `Model.jl` and/or plotter script `Plot.jl` from model 
directory via command-line execution in the background using the function 
`execute_earthbox_script`. This function will run both the model and plotting 
scripts in sequence.

!!! WARNING !!! This function is for Linux/Unix systems only.

# Arguments
- `model_dir::String`
    - Directory from where the model action is taken and that contains Runit.jl, 
       Model.jl and Plot.jl.
- `eb_project_path_from_script::Union{String, Nothing} = nothing`
    - Path to the EarthBox project from the script executing the action. This
       will be overridden by the EarthBox project path specified via command line 
       arguments.
- `root_path_output_from_script::Union{String, Nothing} = nothing`
    - Path to the root output directory from the script executing the action. 
       This will be overridden by the root output directory specified via command 
       line arguments.

# Command Line Arguments
- `case_name=<case_name>`
    - The name of the model case defined in `CaseInputs.jl`. The case name is 
       used when running the model to select inputs from the model case 
       collection defined in `CaseInputs.jl`. The case name is also used to define
       the output path in customizable path functions in `Paths.jl` and `Plot.jl`.
- `root_path=<path/to/root/output>`
    - The path to the root output directory. This will be override the root 
       output directory specified via function argument.
- `eb_path=<path/to/earthbox/project>`
    - The path to the EarthBox project used when running `Model.jl` and `Plot.jl`. 
       If `eb_path` is not defined it will be assumed that EarthBox is 
       installed in a standard Julia location or the `JULIA_PROJECT` environment 
       variable is set equal to the EarthBox project path.
- `model_output_path=<path/to/model/output>`
    - The path to the model output directory. 
- `istart=<istart>`
    - The starting time step for plotting.
- `iend=<iend>`
    - The ending time step for plotting.
- `run_model`
    - Action command for running the model. 
- `plot_markers`
    - Action command for plotting the markers.
- `plot_scalars`
    - Action command for plotting the scalars.

Usage:

Assuming that `run_earthbox()` is called from a script called Runit.jl from
a model directory that contains Model.jl and Plot.jl...

To run model, plot markers, plot scalars and define all paths on command line, 
run:

```bash
julia Runit.jl case_name=<case_name> run_model plot_markers plot_scalars istart=1 iend=10 eb_path=<eb_path> model_output_path=<path/to/model/output> 
```

or to use `root_path` instead of `model_output_path`, run:

```bash
julia Runit.jl case_name=<case_name> run_model plot_markers plot_scalars istart=1 iend=10 root_path=<path/to/root/output>` 
```

To run model, plot markers and plot scalars and use path definitions in Model.jl 
and Plot.jl, run:

```bash
julia Runit.jl case_name=<case_name> run_model plot_markers plot_scalars istart=1 iend=10
```

To run model only and use path definitions in Model.jl and Plot.jl, run:

```
julia Runit.jl case_name=<case_name> run_model
```

If a case name is not specified, the default case name is "case0".

"""
function run_earthbox(;
    model_dir::String,
    eb_project_path_from_script::Union{String, Nothing} = nothing,
    root_path_output_from_script::Union{String, Nothing} = nothing
)::Nothing
    model_case_name = get_case_name_from_cl_args()
    root_path_output = manage_root_path_output(root_path_output_from_script)
    model_output_path = GetPaths.get_model_output_path(model_case_name, root_path_output)
    eb_path = manage_earthbox_project_path(eb_project_path_from_script)
    run_model, plot_markers, plot_scalars = get_runit_actions_from_cl_args()
    check_paths(eb_path, model_dir)
    istart, iend = define_istart_and_iend()
    print_info(
        model_case_name, eb_path, run_model, plot_markers, 
        plot_scalars, istart, iend
    )
    if run_model
        println(">> Running model...")
        execute_earthbox_script(
            model_dir          = model_dir,
            eb_path            = eb_path,
            command_type       = "model",
            model_case_name    = model_case_name,
            model_logfile_name = "model_$(model_case_name).log",
            model_output_path  = model_output_path
        )
    end
    if plot_scalars
        println(">> Plotting scalars...")
        command_type = "scalar_plots"
        check_plotter_command(command_type)
        execute_earthbox_script(
            model_dir         = model_dir,
            eb_path           = eb_path,
            command_type      = command_type,
            model_case_name   = model_case_name,
            istart            = istart,
            iend              = iend,
            model_output_path = model_output_path
        )
    end
    if plot_markers
        println(">> Plotting markers...")
        command_type = "marker_plots"
        check_plotter_command(command_type)
        execute_earthbox_script(
            model_dir         = model_dir,
            eb_path           = eb_path,
            command_type      = command_type,
            model_case_name   = model_case_name,
            istart            = istart,
            iend              = iend,
            model_output_path = model_output_path
        )
    end
    return nothing
end

function manage_root_path_output(
    root_path_output_from_script::Union{String, Nothing}
)::Union{String, Nothing}
    root_path_output = get_root_path_from_args()
    # If the root path from cl is nothing then use the root path from the script
    if isnothing(root_path_output)
        root_path_output = root_path_output_from_script
    end
    return root_path_output
end

function define_istart_and_iend()::Tuple{Int64, Int64}
    istart = get_istart()
    iend = get_iend()
    if iend === nothing
        iend = istart
    end
    return istart, iend
end

function get_model_output_path(
    get_model_output_path_user_func::Union{Function, Nothing} = nothing
)::Union{String, Nothing}
    if !isnothing(get_model_output_path_user_func)
        # Here we use the user defined function for defining the model path to
        # define a default path that is used if no command line path is found
        model_output_path = get_model_output_path_from_cl_args(get_model_output_path_user_func())
    else
        model_output_path = nothing
    end
    return model_output_path
end

function check_paths(
    eb_path::Union{String, Nothing},
    model_dir::String
)::Nothing
    if eb_path !== nothing
        if !isdir(eb_path)
            throw(ArgumentError("EarthBox project path not found: $(eb_path)"))
        end
    end
    if !isdir(model_dir)
        throw(ArgumentError("Model directory not found: $(model_dir)"))
    end
    # Check if Model.jl exists in model_dir
    model_jl_path = joinpath(model_dir, "Model.jl")
    if !isfile(model_jl_path)
        throw(ArgumentError("Model.jl not found in model directory: $(model_jl_path)"))
    end
    # Check if Plot.jl exists in model_dir
    plot_jl_path = joinpath(model_dir, "Plot.jl")
    if !isfile(plot_jl_path)
        throw(ArgumentError("Plot.jl not found in model directory: $(plot_jl_path)"))
    end
    return nothing
end

function print_info(
    model_case_name::String,
    eb_path::Union{String, Nothing},
    run_model::Bool,
    plot_markers::Bool,
    plot_scalars::Bool,
    istart::Int64,
    iend::Int64
)::Nothing
    println(">> Model case name: $model_case_name")
    println(">> EarthBox project path: $eb_path")
    println(">> Run model: $run_model")
    println(">> Plot markers: $plot_markers")
    println(">> Plot scalars: $plot_scalars")
    if plot_markers || plot_scalars
        println(">> istart: $istart")
        println(">> iend: $iend")
    end
    return nothing
end

""" 
    execute_earthbox_script(;
        model_dir::String,
        eb_path::Union{String, Nothing},
        command_type::String,
        model_case_name::String,
        model_logfile_name::String = "model.log",
        model_output_path::Union{String, Nothing} = nothing
        istart::Int64 = 1,
        iend::Int64 = 1,
    )

Execute earthbox scripts via command line.

!!! WARNING !!! This function is for Linux/Unix systems only.

Processes are executed via command line using the Julia `run` function as follows:
- `run(bash -c <command>, wait = true)`

where `<command>` takes the form of:
- `julia [--project=<eb_path>] --startup-file=no <model_dir>/Model.jl case_name=<case_name> model_output_path=<model_output_path> > <model_logfile_name>`
- `julia [--project=<eb_path>] --startup-file=no <model_dir>/Plot.jl marker_plots istart=<istart> iend=<iend> model_case_name=<case_name> model_output_path=<model_output_path> > marker_plots_<model_logfile_name>`
- `julia [--project=<eb_path>] --startup-file=no <model_dir>/Plot.jl scalar_plots istart=<istart> iend=<iend> model_case_name=<case_name> model_output_path=<model_output_path> > scalar_plots_<model_logfile_name>`

If `eb_path` is not specified, then it is assumed that EarthBox is installed in 
a standard Julia location or the `JULIA_PROJECT` environment variable is set equal 
to the EarthBox project path. In this case the `--project=<eb_path>` argument is not 
included in the command.

If `model_output_path` is not specified, then it is assumed that the model output path
is defined internally within the Model.jl and Plot.jl scripts.

The following scripts are expected to be in the model directory:
- `Model.jl`
    - Model script that initializes the model and runs time steps via the EarthBox API.
       This script must be setup to be executed via command line with the command line 
       arguments `case_name=<case_name>`.
- `Plot.jl`
    - Plot script that plots the model results using the EarthBox API. This script 
       must be setup to be executed via command line with the command line arguments 
       `istart=<istart> iend=<iend> model_case_name=<case_name>`.

# Arguments
- `model_dir::String`
    - Directory from where the script action is taken and that contains the script.
- `eb_path::Union{String, Nothing}`
    - EarthBox project path
- `command_type::String`
    - Command type
        - "model": Run the model script `Model.jl`
        - "marker_plots": Plot markers using the `Plot` module from a Plot.jl script
        - "scalar_plots": Plot scalars using the `Plot` module from a Plot.jl script
- `model_case_name::String`
    - Model case name like "case0", "case1", "case2", etc.
- `model_logfile_name::String`
    - Model log file name (only used if command_type is "model")
- `model_output_path::Union{String, Nothing}`
    - Model output path where model output will be saved and log files will be copied
- `istart::Int64`
    - Starting time step. Default is 1.
- `iend::Int64`
    - Ending time step. Default is 1.
"""
function execute_earthbox_script(;
    model_dir::String,
    eb_path::Union{String, Nothing},
    command_type::String,
    model_case_name::String,
    model_logfile_name::String = "model.log",
    model_output_path::Union{String, Nothing} = nothing,
    istart::Int64 = 1,
    iend::Int64 = 1
)::Nothing
    command = build_command(
        model_dir, eb_path, command_type, model_case_name, 
        model_logfile_name, model_output_path, istart, iend
        )
    try
        println(">> Executing: $(command)")
        result = run(`bash -c "$(command)"`, wait = true)
        if result.exitcode != 0
            println(">> Command failed with exit code: $(result.exitcode)")
            return
        end
        println(">> Command executed successfully")
    catch e
        println(">> !!! ERROR !!! An error occurred while executing '$(command)': $(e)")
        return
    end
    if command_type == "model" && model_output_path !== nothing
        copy_output_log(model_logfile_name, model_output_path)
    end
    return nothing
end

function build_command(
    model_dir::String,
    eb_path::Union{String, Nothing},
    command_type::String,
    model_case_name::String,
    model_logfile_name::String,
    model_output_path::Union{String, Nothing},
    istart::Int64,
    iend::Int64
)::String
    command_registry = get_command_registry(model_dir, eb_path)
    base_command = get(command_registry, command_type, nothing)
    security_checks(base_command, command_registry, model_logfile_name)
    
    model_output_path_part = get_model_output_path_part(model_output_path)
    case_name_part = get_case_name_part(model_case_name)
    istart_iend_part = get_istart_iend_part(istart, iend)
    
    if command_type == "model"
        command = "$(base_command) $(case_name_part) $(model_output_path_part) > $(model_logfile_name)"
    elseif command_type == "marker_plots"
        command = "$(base_command) $(istart_iend_part) $(case_name_part) $(model_output_path_part) > marker_plots_$(model_logfile_name)"
    elseif command_type == "scalar_plots"
        command = "$(base_command) $(istart_iend_part) $(case_name_part) $(model_output_path_part) > scalar_plots_$(model_logfile_name)"
    end
    return command
end

function get_command_registry(
    model_dir::String, 
    eb_path::Union{String, Nothing}
)::Dict{String, String}
    if eb_path === nothing
        # No project path specified, assume EarthBox is installed in a standard Julia location or 
        # the JULIA_PROJECT environment variable is set equal to the EarthBox project path.
        eb_part = ""
    else
        eb_part = " --project=$(eb_path)"
    end
    command_registry = Dict{String, String}(
        "model"        => "julia$(eb_part) --startup-file=no $(model_dir)/Model.jl",
        "marker_plots" => "julia $(eb_part) --startup-file=no $(model_dir)/Plot.jl marker_plots",
        "scalar_plots" => "julia $(eb_part) --startup-file=no $(model_dir)/Plot.jl scalar_plots",
    )
    return command_registry
end

function security_checks(
    command::Union{String, Nothing},
    command_registry::Dict{String, String},
    model_logfile_name::String
)::Nothing
    check_command(command, command_registry)
    
    if !isa(model_logfile_name, String)
        throw(ArgumentError("Model log file name must be a string."))
    end
    
    return nothing
end

function check_command(
    command::Union{String, Nothing},
    command_registry::Dict{String, String}
)::Nothing
    if command === nothing
        println("Registered commands:")
        for (key, cmd) in command_registry
            println("$(key) : $(cmd)")
        end
        throw(ArgumentError("Command '$(command)' not found in command registry."))
    end
    
    return nothing
end

function get_case_name_part(model_case_name::String)::String
    return "case_name=$(model_case_name)"
end

function get_model_output_path_part(model_output_path::Union{String, Nothing})::String
    if model_output_path !== nothing
        return "model_output_path=$(model_output_path)"
    else
        return ""
    end
end

function get_istart_iend_part(istart::Int64, iend::Int64)::String
    return "istart=$(istart) iend=$(iend)"
end

function check_plotter_command(plotter_command::String)::Nothing
    valid_commands = ["marker_plots", "scalar_plots"]
    if !(plotter_command in valid_commands)
        throw(ArgumentError("Invalid plotter command: $(plotter_command)"))
    end
end


function copy_output_log(
    output_file_name::String,
    final_log_dir_path::String
)::Nothing
    if isfile(output_file_name)
        try
            mv(output_file_name, joinpath(final_log_dir_path, output_file_name))
            println(">> Moved $(output_file_name) to $(joinpath(final_log_dir_path, output_file_name))")
        catch e
            println(">> !!! WARNING !!! An error occurred while moving $(output_file_name) to $(final_log_dir_path): $(e)")
        end
    else
        println(">> $(output_file_name) not found. Ensure your commands generate this file.")
    end
    
    return nothing
end

""" 
    remote_model_loop(;
        models::Dict{String, Vector{String}},
        model_dir_path::String,
        run_model::Union{Bool, Nothing} = nothing,
        plot_markers::Union{Bool, Nothing} = nothing,
        plot_scalars::Union{Bool, Nothing} = nothing,
        istart::Union{Int64, Nothing} = nothing,
        iend::Union{Int64, Nothing} = nothing
    )

Execute a script on a remote machine in the background by looping over the
function `execute_remote_script_in_background` for each model case.

!!! WARNING !!! This function is for Linux/Unix systems only.

# Arguments
- `models::Dict{String, Vector{String}}`
    - A dictionary of host names and model case names like:
        - "hostname1" => ["case0", "case1", "case2"]
        - "hostname2" => ["case3", "case4", "case5"]
- `model_dir_path::String`
    - The path to the directory that contains the model directories
- `run_model::Union{Bool, Nothing}`
    - Run the model (default: nothing)
- `plot_markers::Union{Bool, Nothing}`
    - Plot markers (default: nothing)
- `plot_scalars::Union{Bool, Nothing}`
    - Plot scalars (default: nothing)
- `istart::Union{Int64, Nothing}`
    - Start index for plotting (default: nothing)
- `iend::Union{Int64, Nothing}`
    - End index for plotting (default: nothing)
"""
function remote_model_loop(
    models::Dict{String, Vector{String}},
    model_dir_path::String;
    run_model::Union{Bool, Nothing} = nothing,
    plot_markers::Union{Bool, Nothing} = nothing,
    plot_scalars::Union{Bool, Nothing} = nothing,
    istart::Union{Int64, Nothing} = nothing,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    for (hostname, model_case_names) in models
        println(">> On hostname $(hostname), executing scripts for the following model cases: $(model_case_names)")
        for model_case_name in model_case_names
            print_execution_msg(hostname, model_case_name)
            execute_remote_script_in_background(
                hostname,
                model_dir_path,
                model_case_name,
                run_model,
                plot_markers,
                plot_scalars,
                istart,
                iend
            )
            # Wait for a few seconds before starting the next model
            sleep(2)
        end
    end
    
    return nothing
end

""" 
    local_model_loop(;
        models::Dict{String, Vector{String}},
        model_dir_path::String,
        run_model::Union{Bool, Nothing} = nothing,
        plot_markers::Union{Bool, Nothing} = nothing,
        plot_scalars::Union{Bool, Nothing} = nothing,
        istart::Union{Int64, Nothing} = nothing,
        iend::Union{Int64, Nothing} = nothing
    )

Execute a script locally by looping over the
function `execute_remote_script_in_background` for each model case.
    
!!! WARNING !!! This function is for Linux/Unix systems only.

# Arguments
- `models::Dict{String, Vector{String}}`
    - A dictionary of model directories containing `Model.jl` and `Plot.jl` scripts 
        and model case names like:
        - "/path/to/model_dir1" => ["case0", "case1", "case2"]
        - "/path/to/model_dir2" => ["case3", "case4", "case5"]
- `run_model::Union{Bool, Nothing}`
    - Run the model (default: nothing)
- `plot_markers::Union{Bool, Nothing}`
    - Plot markers (default: nothing)
- `plot_scalars::Union{Bool, Nothing}`
    - Plot scalars (default: nothing)
- `istart::Union{Int64, Nothing}`
    - Start index for plotting (default: nothing)
- `iend::Union{Int64, Nothing}`
    - End index for plotting (default: nothing)
"""
function local_model_loop(;
    models::Dict{String, Vector{String}},
    run_model::Union{Bool, Nothing} = nothing,
    plot_markers::Union{Bool, Nothing} = nothing,
    plot_scalars::Union{Bool, Nothing} = nothing,
    istart::Union{Int64, Nothing} = nothing,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    for (model_dir, model_case_names) in models
        println(">> In model directory $(model_dir), executing scripts for the following model cases: $(model_case_names)")
        # Wait for a few seconds before starting the next model
        for model_case_name in model_case_names
            println("    >> Working on model case: $(model_case_name)")
            execute_local_script_in_background(
                model_dir,
                model_case_name,
                run_model,
                plot_markers,
                plot_scalars,
                istart,
                iend
            )
            sleep(2)
        end
    end
    
    return nothing
end

function print_execution_msg(
        hostname::String,
        model_case_name::String
)::Nothing
    println(">> On hostname $(hostname), executing script for $(model_case_name) in background.")
    return nothing
end

function print_execution_msg_local(
        model_dir::String,
        model_case_name::String
)::Nothing
    println(">> In model directory $(model_dir), executing script for $(model_case_name).")
    return nothing
end

""" 
    execute_remote_script_in_background(;
        hostname::String,
        model_dir_path::String,
        model_case_name::String,
        run_model::Union{Bool, Nothing} = nothing,
        plot_markers::Union{Bool, Nothing} = nothing,
        plot_scalars::Union{Bool, Nothing} = nothing,
        istart::Union{Int64, Nothing} = nothing,
        iend::Union{Int64, Nothing} = nothing
    )

This function runs the `Runit.jl` script on a remote machine in the background as follows:

    ssh <hostname> 'cd <model_dir_path> && nohup julia --startup=no <model_dir_path>/Runit.jl case_name=<case_name> <runit_actions> <istart_iend> > /dev/null 2>&1 &'

where:
- `<hostname>` is the name of the remote machine.
- `<model_dir_path>` is the path to the directory that contains the model input files and the 
    `Runit.jl` script.
- `<model_case_name>` is the name of the model case like "case0", "case1", "case2", etc.
- `<runit_actions>` are the actions to be taken by the `Runit.jl` script like "run_model", 
     "plot_markers", "plot_scalars" or a combination of these separated by spaces.
- `<istart_iend>` are the starting and ending time steps. For example, "istart=1 iend=10".

The `Runit.jl` script should be in a model directory that contains `Model.jl` 
and `Plot.jl` scripts that are properly configured and setup for command-line 
execution of the `run_earthbox()` function. For example, the `Runit.jl` script 
should contain the following code or similar code that achieves the same result:

```julia
module Runit

using EarthBox
include("Model.jl")
import .Model: EB_PROJECT_PATH, ROOT_PATH_OUTPUT
run_earthbox(
    model_dir = PATH_TO_MODEL_DIRECTORY,
    eb_project_path_from_script = EB_PROJECT_PATH,
    root_path_output_from_script = ROOT_PATH_OUTPUT
    )

end

if abspath(PROGRAM_FILE) == @__FILE__
    Runit.main()
end
```

# Arguments
- `hostname::String`
    - The name of the remote machine.
- `model_dir_path::String`
    - The path to the directory that contains the model input files and the 
       `Runit.jl` script.
- `model_case_name::String`
    - The name of the model case. This is the name of the case that will be 
       executed as defined in the model script like "case0", "case1", "case2", etc.
- `run_model::Union{Bool, Nothing}`
    - Run the model (default: nothing)
- `plot_markers::Union{Bool, Nothing}`
    - Plot markers (default: nothing)
- `plot_scalars::Union{Bool, Nothing}`
    - Plot scalars (default: nothing)
- `istart::Union{Int64, Nothing}`
    - Starting time step. Default is 1.
- `iend::Union{Int64, Nothing}`
    - Ending time step. Default is nothing.

It is assumed that passwordless SSH access is configured for the remote machine and
hostnames are resolved to IP addresses.

# Setting up passwordless SSH access

On Ubuntu/PopOS/Debian-based systems you can map IP addresses to hostnames in 
the `/etc/hosts` file by editing the file with the following command:

```bash
sudo vi /etc/hosts
```

and adding an IP-address and hostname pair on each line for each remote machine.

Passwordless SSH access can be configured on Ubuntu/PopOS/Debian-based systems by 
following the instructions:

[1] Generate SSH key files:

```bash
ssh-keygen -t rsa
```

Use `id_rsa` for the private key that will be located at `/home/username/.ssh` 
and follow the prompts to enter a passphrase or use an empty passphrase. A
public key file will also be created that will be shared with the remote 
machine (see step [2] below).

[2] Copy the public SSH key to the remote machine:

```bash
ssh-copy-id <username>@<hostname>
```

"""
function execute_remote_script_in_background(
    hostname::String,
    model_dir_path::String,
    model_case_name::String,
    run_model::Union{Bool, Nothing} = nothing,
    plot_markers::Union{Bool, Nothing} = nothing,
    plot_scalars::Union{Bool, Nothing} = nothing,
    istart::Union{Int64, Nothing} = nothing,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    runit_actions = get_runit_actions(run_model, plot_markers, plot_scalars)
    script_name = "Runit.jl"
    script_path = joinpath(model_dir_path, script_name)
    try
        ssh_command = "ssh $(hostname) 'cd $(model_dir_path) && nohup julia --startup=no $(script_path) case_name=$(model_case_name) $(runit_actions) $(get_istart_iend(istart, iend)) > /dev/null 2>&1 &'"
        println(">> Executing: $(ssh_command)")
        run(`bash -c "$(ssh_command)"`, wait = false)
    catch e
        println(">> !!! ERROR !!! An error occurred: $(e)")
    end
    
    return nothing
end

"""
    execute_local_script_in_background(;
        model_dir_path::String,
        model_case_name::String,
        run_model::Union{Bool, Nothing} = nothing,
        plot_markers::Union{Bool, Nothing} = nothing,
        plot_scalars::Union{Bool, Nothing} = nothing,
        istart::Union{Int64, Nothing} = nothing,
        iend::Union{Int64, Nothing} = nothing
    )

Execute the `Runit.jl` script with a call to `run_earthbox()` locally in the 
background. The `Runit.jl` script should be in a model directory that contains 
`Model.jl` and `Plot.jl` scripts that are properly configured and setup for 
command-line execution of the `run_earthbox()` function. For example, the 
`Runit.jl` script should contain the following code or similar code that 
achieves the same result:

```julia
module Runit

using EarthBox
include("Model.jl")
import .Model: EB_PROJECT_PATH, ROOT_PATH_OUTPUT
run_earthbox(
    model_dir = PATH_TO_MODEL_DIRECTORY,
    eb_project_path_from_script = EB_PROJECT_PATH,
    root_path_output_from_script = ROOT_PATH_OUTPUT
    )

end

if abspath(PROGRAM_FILE) == @__FILE__
    Runit.main()
end
```

This function runs the `Runit.jl` script locally in the background as follows:

    julia --startup=no <model_dir_path>/Runit.jl case_name=<case_name> <runit_actions> <istart_iend> > /dev/null 2>&1 &'

# Arguments
- `model_dir_path::String`
    - The path to the directory that contains the `Model.jl` and `Plot.jl` scripts.
- `model_case_name::String`
    - The name of the model case. This is the name of the case that will be 
       executed as defined in the model script like "case0", "case1", "case2", etc.
- `run_model::Union{Bool, Nothing}`
    - Run the model (default: nothing)
- `plot_markers::Union{Bool, Nothing}`
    - Plot markers (default: nothing)
- `plot_scalars::Union{Bool, Nothing}`
    - Plot scalars (default: nothing)
- `istart::Union{Int64, Nothing}`
    - Starting time step. Default is 1.
- `iend::Union{Int64, Nothing}`
    - Ending time step. Default is nothing.

"""
function execute_local_script_in_background(
    model_dir_path::String,
    model_case_name::String,
    run_model::Union{Bool, Nothing} = nothing,
    plot_markers::Union{Bool, Nothing} = nothing,
    plot_scalars::Union{Bool, Nothing} = nothing,
    istart::Union{Int64, Nothing} = nothing,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    runit_actions = get_runit_actions(run_model, plot_markers, plot_scalars)
    script_name = "Runit.jl"
    script_path = joinpath(model_dir_path, script_name)
    try
        local_command = "cd $(model_dir_path) &&  nohup julia --startup=no $(script_path) case_name=$(model_case_name) $(runit_actions) $(get_istart_iend(istart, iend)) > /dev/null 2>&1 &"
        println(">> Executing: $(local_command)")
        run(`bash -c "$(local_command)"`, wait = false)
    catch e
        println(">> !!! ERROR !!! An error occurred: $(e)")
    end
    return nothing
end

function get_runit_actions(
    run_model::Union{Bool, Nothing},
    plot_markers::Union{Bool, Nothing},
    plot_scalars::Union{Bool, Nothing},
)::String
    runit_actions = ""
    if run_model === true
        runit_actions = "run_model"
    end
    if plot_markers === true
        runit_actions = runit_actions * " " * "plot_markers"
    end
    if plot_scalars === true
        runit_actions = runit_actions * " " * "plot_scalars"
    end
    return runit_actions
end

function get_istart_iend(
    istart::Union{Int64, Nothing},
    iend::Union{Int64, Nothing}
)::String
    istart_iend = ""
    if istart !== nothing
        istart_iend = "istart=$(istart)"
    end
    if iend !== nothing
        istart_iend = "$(istart_iend) iend=$(iend)"
    end
    return istart_iend
end

end # module
