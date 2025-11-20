""" Functions used to move files to longer term storage """
module FileMover

import Base.Threads
import Base.Filesystem
import Base.catch_backtrace
using Printf

const PRINT_MOVE_FILE_MESSAGE = true

""" A struct to move output files from fast storage to longer-term storage.

This struct will monitor a directory for files that are no longer being
written to and move them to a destination directory for longer term storage.

# Fields

- `output_directory::Union{String, Nothing}`:
    - The output directory to monitor for files.
- `storage_directory::Union{String, Nothing}`:
    - The directory to move stable files to for longer term storage.
- `check_interval::Union{Int, Float64}`:
    - How often (in seconds) to check for stable files.
- `stable_time::Union{Int, Float64}`:
    - Time (in seconds) a file must remain unchanged to be considered stable.
- `file_mover_task::Union{Task, Nothing}`: 
    - The task that monitors the output directory and moves files.
- `stop_event::Union{Threads.Event, Nothing}`: 
    - An event to signal the file mover task to stop.
- `do_not_move::Vector{String}`: 
    - List of files that should not be moved.
"""
mutable struct FileMoverState
    output_directory::Union{String, Nothing}
    storage_directory::Union{String, Nothing}
    check_interval::Union{Int, Float64}
    stable_time::Union{Int, Float64}
    file_mover_task::Union{Task, Nothing}
    stop_event::Union{Threads.Event, Nothing}
    do_not_move::Vector{String}
end

function FileMoverState(;
    output_directory::Union{String, Nothing}=nothing,
    storage_directory::Union{String, Nothing}=nothing,
    check_interval::Union{Int, Float64}=5,
    stable_time::Union{Int, Float64}=25
)::FileMoverState
    do_not_move = get_list_of_files_to_not_move()
    if storage_directory === nothing || output_directory === nothing
        file_mover_task, stop_event = nothing, nothing
    else
        println(">> Initializing file mover process...")
        file_mover_task, stop_event = initialize_file_mover_process(
            output_directory, storage_directory, check_interval, 
            stable_time, do_not_move
        )
    end
    return FileMoverState(
        output_directory, 
        storage_directory, 
        check_interval, 
        stable_time,
        file_mover_task, 
        stop_event, 
        do_not_move
    )
end

function file_mover_cleanup(file_mover::FileMoverState)::Nothing
    if file_mover.file_mover_task !== nothing && file_mover.stop_event !== nothing
        delay_until_all_files_moved(file_mover)
        stop_file_mover_process(file_mover)
        move_remaining_files(file_mover)
    end
    return nothing
end

""" Get a list of files that should not be moved.

Large backup files that are overwritten during model runs should not be
moved since (1) backup files should be available in fast storage and (2)
these files may be overwritten while being moved to long-term storage
which can cause crashes.
"""
function get_list_of_files_to_not_move()::Vector{String}
    return ["model_backup.jld"]
end

function initialize_file_mover_process(
    output_directory::String,
    storage_directory::String,
    check_interval::Union{Int, Float64},
    stable_time::Union{Int, Float64},
    do_not_move::Vector{String}
)::Tuple{Task, Threads.Event}
    # Create an event to signal the process to stop
    stop_event = Threads.Event()
    # Start the file mover in a separate task
    file_mover_task = Threads.@async move_stable_files(
        output_directory, storage_directory, check_interval, stable_time,
        do_not_move, stop_event
    )
    return file_mover_task, stop_event
end

function stop_file_mover_process(file_mover::FileMoverState)::Nothing
    if file_mover.file_mover_task !== nothing
        # Signal the task to stop
        notify(file_mover.stop_event)
        # Wait for the task to complete
        wait(file_mover.file_mover_task)
    end
    return nothing
end

function delay_until_all_files_moved(file_mover::FileMoverState)::Nothing
    if file_mover.output_directory === nothing || file_mover.storage_directory === nothing
        return nothing
    end
    max_delay = 60.0
    t_start = time()
    dt = 0.0
    file_list = get_list_of_earthbox_output_files(file_mover)
    nfiles = length(file_list)
    while nfiles > 0 && dt < max_delay
        sleep(5)
        file_list = get_list_of_earthbox_output_files(file_mover)
        nfiles = length(file_list)
        t_current = time()
        dt = t_current - t_start
    end
    if dt >= max_delay
        println("***   Delay exceeded 60 seconds. Exiting file mover..")
    else
        println("***   All files have been moved to long-term storage...")
    end
    return nothing
end

function get_list_of_earthbox_output_files(file_mover::FileMoverState)::Vector{String}
    file_list = String[]
    for file_name in Filesystem.readdir(file_mover.output_directory)
        file_path = Filesystem.joinpath(file_mover.output_directory, file_name)
        if Filesystem.isfile(file_path)
            if !(file_name in file_mover.do_not_move)
                push!(file_list, file_name)
            end
        end
    end
    return file_list
end

""" Monitors a directory and moves files that are no longer being written to.
"""
function move_stable_files(
    output_directory::String,
    storage_directory::String,
    check_interval::Union{Int, Float64},
    stable_time::Union{Int, Float64},
    do_not_move::Vector{String},
    stop_event::Threads.Event
)::Nothing
    if !Filesystem.isdir(storage_directory)
        Filesystem.mkpath(storage_directory)
    end
    file_last_modified = Dict{String, Float64}()
    delta_file_sizes = Dict{String, Float64}()
    # Use a different approach to check for stop signal
    # Create a task that waits for the stop event and sets a flag
    stop_flag = Ref(false)
    stop_task = @async begin
        wait(stop_event)
        stop_flag[] = true
    end
    while !stop_flag[]
        current_time = time()
        file_sizes_before_sleep = get_size_of_all_files(output_directory)
        move_all_files_in_directory(
            output_directory, storage_directory, current_time, stable_time,
            file_last_modified, delta_file_sizes, do_not_move
        )
        file_last_modified = check_last_modified_for_file_existence(
            file_last_modified, output_directory
        )
        
        if !stop_flag[]
            sleep(check_interval)
            file_sizes_after_sleep = get_size_of_all_files(output_directory)
            delta_file_sizes = calculate_delta_files_sizes(
                file_sizes_before_sleep, file_sizes_after_sleep
            )
        end
    end
    return nothing
end

""" Calculate the change in file sizes between two time points.
"""
function calculate_delta_files_sizes(
    file_sizes_before_sleep::Dict{String, Float64},
    file_sizes_after_sleep::Dict{String, Float64}
)::Dict{String, Float64}
    delta_file_sizes = Dict{String, Float64}()
    for file_name in keys(file_sizes_before_sleep)
        if haskey(file_sizes_after_sleep, file_name)
            delta_file_sizes[file_name] = (
                file_sizes_after_sleep[file_name] - file_sizes_before_sleep[file_name]
            )
        end
    end
    return delta_file_sizes
end

""" Get the size of all files in a directory.
"""
function get_size_of_all_files(output_directory::String)::Dict{String, Float64}
    file_sizes = Dict{String, Float64}()
    for file_name in Filesystem.readdir(output_directory)
        file_path = Filesystem.joinpath(output_directory, file_name)
        if Filesystem.isfile(file_path)
            file_size = Filesystem.filesize(file_path)
            file_sizes[file_name] = file_size
        end
    end
    return file_sizes
end

""" Move files that are no longer being written to.
"""
function move_all_files_in_directory(
    output_directory::String,
    storage_directory::String,
    current_time::Float64,
    stable_time::Union{Int, Float64},
    file_last_modified::Dict{String, Float64},
    delta_file_sizes::Dict{String, Float64},
    do_not_move::Vector{String}
)::Nothing
    for file_name in Filesystem.readdir(output_directory)
        file_path = Filesystem.joinpath(output_directory, file_name)
        if Filesystem.isfile(file_path)
            last_modified = Filesystem.mtime(file_path)
            delta_size = get(delta_file_sizes, file_name, nothing)
            if haskey(file_last_modified, file_name) && !(file_name in do_not_move)
                delta_time = current_time - last_modified
                check_and_move_file(
                    file_path, file_name, storage_directory, delta_time,
                    stable_time, delta_size, file_last_modified
                )
            else
                file_last_modified[file_name] = last_modified
            end
        end
    end
    return nothing
end

function check_and_move_file(
    file_path::String,
    file_name::String,
    storage_directory::String,
    delta_time::Float64,
    stable_time::Union{Int, Float64},
    delta_size::Union{Float64, Nothing},
    file_last_modified::Dict{String, Float64}
)::Nothing
    # If delta size has not been measured yet do not move file
    if delta_size === nothing
        return nothing
    end
    # Only move file if file has not been modified for stable_time seconds
    # and the file size has not changed since the last check
    if delta_time >= stable_time && delta_size == 0
        # Re-check if file still exists before moving to avoid race conditions
        if Filesystem.isfile(file_path)
            try
                Filesystem.mv(file_path, Filesystem.joinpath(storage_directory, file_name); force=true)
                print_move_file_message(file_name, storage_directory)
                # Remove file from last modified dictionary
                delete!(file_last_modified, file_name)
            catch e
                handle_exception(e)
            end
        else
            # Remove file from last modified dictionary
            delete!(file_last_modified, file_name)
        end
    end
    return nothing
end

""" Create a fresh last modified dict only with files that exist.

Only files that still exist in the source directory are included in the
updated dictionary.
"""
function check_last_modified_for_file_existence(
    file_last_modified_old::Dict{String, Float64},
    output_directory::String
)::Dict{String, Float64}
    file_last_modified = Dict{String, Float64}()
    for (file_name, last_modified) in file_last_modified_old
        if Filesystem.isfile(Filesystem.joinpath(output_directory, file_name))
            file_last_modified[file_name] = last_modified
        end
    end
    return file_last_modified
end

""" Handle an exception.
"""
function handle_exception(exception::Exception)::Nothing
    println("An unexpected error occurred in file mover:")
    println("Type: ", typeof(exception))
    println("Message: ", exception)
    println("Detailed Traceback:")
    showerror(stdout, exception, catch_backtrace())
    return nothing
end

function move_remaining_files(file_mover::FileMoverState)::Nothing
    if file_mover.output_directory === nothing || file_mover.storage_directory === nothing
        return nothing
    end
    for file_name in Filesystem.readdir(file_mover.output_directory)
        file_path = Filesystem.joinpath(file_mover.output_directory, file_name)
        if Filesystem.isfile(file_path)
            if !(file_name in file_mover.do_not_move)
                move_single_file(file_path, file_name, file_mover.storage_directory)
            end
        end
    end
    return nothing
end

function move_single_file(
    file_path::String,
    file_name::String,
    storage_directory::String
)::Nothing
    # Check if file still exists before moving to avoid race conditions
    if Filesystem.isfile(file_path)
        try
            Filesystem.mv(file_path, Filesystem.joinpath(storage_directory, file_name); force=true)
            print_move_file_message(file_name, storage_directory)
        catch e
            handle_exception(e)
        end
    end
    return nothing
end

function print_move_file_message(file_name::String, storage_directory::String)::Nothing
    if PRINT_MOVE_FILE_MESSAGE
        println("   ---> Moved file ", file_name, " to ", storage_directory)
    end
    return nothing
end

end # module 