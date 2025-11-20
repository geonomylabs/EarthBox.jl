module Logging

using Printf
using Dates

# Global debug flag
const DEBUG = false

# Log levels
@enum LogLevel begin
    DEBUG = 1
    INFO = 2
    WARN = 3
    ERROR = 4
    NONE = 5
end

# Global state
mutable struct LoggingState
    current_level::LogLevel
    output_stream::IO
    prefix::String
    timestamp::Bool
end

# Default state
const DEFAULT_STATE = LoggingState(INFO, stdout, "EarthBox", true)

# Global state instance
const GLOBAL_STATE = deepcopy(DEFAULT_STATE)

"""
    set_log_level!(level::LogLevel)

Set the global log level. Messages with a level lower than this will not be printed.
"""
function set_log_level!(level::LogLevel)
    GLOBAL_STATE.current_level = level
end

"""
    set_output_stream!(stream::IO)

Set the output stream for logging messages.
"""
function set_output_stream!(stream::IO)
    GLOBAL_STATE.output_stream = stream
end

"""
    set_log_prefix!(prefix::String)

Set the prefix for log messages.
"""
function set_log_prefix!(prefix::String)
    GLOBAL_STATE.prefix = prefix
end

"""
    enable_timestamps!()

Enable timestamps in log messages.
"""
function enable_timestamps!()
    GLOBAL_STATE.timestamp = true
end

"""
    disable_timestamps!()

Disable timestamps in log messages.
"""
function disable_timestamps!()
    GLOBAL_STATE.timestamp = false
end

"""
    reset_logging!()

Reset logging configuration to defaults.
"""
function reset_logging!()
    global GLOBAL_STATE
    GLOBAL_STATE = deepcopy(DEFAULT_STATE)
end

"""
    log(level::LogLevel, message::String; kwargs...)

Log a message at the specified level. Additional keyword arguments will be printed after the message.
"""
function log(level::LogLevel, message::String; kwargs...)
    if level < GLOBAL_STATE.current_level
        return
    end

    timestamp = GLOBAL_STATE.timestamp ? string("[", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "] ") : ""
    prefix = GLOBAL_STATE.prefix
    level_str = string(level)
    
    # Print the main message
    print(GLOBAL_STATE.output_stream, "$timestamp$prefix [$level_str] $message")
    
    # Print any additional keyword arguments
    if !isempty(kwargs)
        print(GLOBAL_STATE.output_stream, " (")
        for (i, (k, v)) in enumerate(kwargs)
            print(GLOBAL_STATE.output_stream, "$k=$v")
            if i < length(kwargs)
                print(GLOBAL_STATE.output_stream, ", ")
            end
        end
        print(GLOBAL_STATE.output_stream, ")")
    end
    
    println(GLOBAL_STATE.output_stream)
    flush(GLOBAL_STATE.output_stream)
end

"""
    debug(message::String; kwargs...)

Log a debug message.
"""
function debug(message::String; kwargs...)
    log(DEBUG, message; kwargs...)
end

"""
    info(message::String; kwargs...)

Log an info message.
"""
function info(message::String; kwargs...)
    log(INFO, message; kwargs...)
end

"""
    warn(message::String; kwargs...)

Log a warning message.
"""
function warn(message::String; kwargs...)
    log(WARN, message; kwargs...)
end

"""
    error(message::String; kwargs...)

Log an error message.
"""
function error(message::String; kwargs...)
    log(ERROR, message; kwargs...)
end

"""
    debug_print(message::String; kwargs...)

Print a debug message if DEBUG is true. This is a simpler alternative to the full logging system
for controlling debug output inside functions.

# Arguments
- `message::String`: The debug message to print
- `kwargs...`: Optional key-value pairs to include in the output

# Example
```julia
debug_print("Processing data", size=size(data), type=typeof(data))
```
"""
function debug_print(message::String; kwargs...)
    if !DEBUG
        return
    end
    
    # Print the main message
    print("DEBUG: $message")
    
    # Print any additional keyword arguments
    if !isempty(kwargs)
        print(" (")
        for (i, (k, v)) in enumerate(kwargs)
            print("$k=$v")
            if i < length(kwargs)
                print(", ")
            end
        end
        print(")")
    end
    
    println()
    flush(stdout)
end

end # module 