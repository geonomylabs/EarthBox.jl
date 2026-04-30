module ArrayMacros

export @arrays

"""
    @arrays (key1 = ArrayData(units, type, grid_type, description),
             key2 = ArrayData(units, type, grid_type, description),
             ...)

Build a NamedTuple of `ArrayData` entries without writing each array's name twice.

Each entry's left-hand-side symbol is automatically spliced into the constructor
call as the `name` argument (position 1), so the constructor calls inside the
macro take only `(units, type, grid_type, description)` instead of
`(name, units, type, grid_type, description)`.

Only `ArrayData` is accepted as a constructor. Each call must pass exactly 4
positional arguments; a 5-argument call is rejected because it indicates the
caller passed `name` explicitly, defeating the purpose of the macro.
"""
macro arrays(ex)
    if !(ex isa Expr) || ex.head !== :tuple
        error("@arrays expects a NamedTuple-style tuple expression, e.g. " *
              "`@arrays (a = ArrayData(\"m\", T, \"basic\", \"desc\"),)`. " *
              "Got: $(ex)")
    end

    rewritten_entries = Any[]
    for entry in ex.args
        if !(entry isa Expr) || entry.head !== :(=)
            error("@arrays: each entry must have the form " *
                  "`name = ArrayData(units, type, grid_type, description)`. " *
                  "Got: $(entry)")
        end

        key = entry.args[1]
        call = entry.args[2]

        if !(key isa Symbol)
            error("@arrays: each entry's key must be a plain symbol. Got: $(key)")
        end

        if !(call isa Expr) || call.head !== :call
            error("@arrays: right-hand side of `$(key)` must be a constructor " *
                  "call. Got: $(call)")
        end

        ctor = call.args[1]
        if !(ctor isa Symbol) || ctor !== :ArrayData
            error("@arrays: only ArrayData constructors are allowed. " *
                  "Got `$(ctor)` for entry `$(key)`.")
        end

        ctor_args = call.args[2:end]
        if length(ctor_args) != 4
            error("@arrays: `ArrayData` for entry `$(key)` must take exactly 4 " *
                  "arguments (units, type, grid_type, description); the macro " *
                  "inserts the name automatically. Got $(length(ctor_args)) " *
                  "arguments. If you wrote the name explicitly as the first " *
                  "argument, remove it.")
        end

        name_str = String(key)
        new_call = Expr(:call, ctor, name_str,
                        ctor_args[1], ctor_args[2], ctor_args[3], ctor_args[4])
        push!(rewritten_entries, Expr(:(=), key, new_call))
    end

    return esc(Expr(:tuple, rewritten_entries...))
end

end # module
