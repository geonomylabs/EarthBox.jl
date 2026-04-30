module ParamMacros

export @params

"""
    @params (key1 = ParameterXxx(value, units, description),
             key2 = ParameterXxx(value, units, description),
             ...)

Build a NamedTuple of parameter objects without writing each parameter's name twice.

Each entry's left-hand-side symbol is automatically spliced into the constructor
call as the `name` argument (position 2), so the constructor calls inside the
macro take only `(value, units, description)` instead of
`(value, name, units, description)`.

Only `ParameterFloat`, `ParameterInt`, and `ParameterStr` are accepted as
constructors. Each call must pass exactly 3 positional arguments
(`value`, `units`, `description`); a 4-argument call is rejected because it
indicates the caller passed `name` explicitly, defeating the purpose of the macro.
"""
macro params(ex)
    if !(ex isa Expr) || ex.head !== :tuple
        error("@params expects a NamedTuple-style tuple expression, e.g. " *
              "`@params (a = ParameterInt(0, \"None\", \"desc\"),)`. " *
              "Got: $(ex)")
    end

    rewritten_entries = Any[]
    for entry in ex.args
        if !(entry isa Expr) || entry.head !== :(=)
            error("@params: each entry must have the form " *
                  "`name = ParameterXxx(value, units, description)`. " *
                  "Got: $(entry)")
        end

        key = entry.args[1]
        call = entry.args[2]

        if !(key isa Symbol)
            error("@params: each entry's key must be a plain symbol. Got: $(key)")
        end

        if !(call isa Expr) || call.head !== :call
            error("@params: right-hand side of `$(key)` must be a constructor " *
                  "call. Got: $(call)")
        end

        ctor = call.args[1]
        if !(ctor isa Symbol) || !(ctor in (:ParameterFloat, :ParameterInt, :ParameterStr))
            error("@params: only ParameterFloat / ParameterInt / ParameterStr " *
                  "constructors are allowed. Got `$(ctor)` for entry `$(key)`.")
        end

        ctor_args = call.args[2:end]
        if length(ctor_args) != 3
            error("@params: `$(ctor)` for entry `$(key)` must take exactly 3 " *
                  "arguments (value, units, description); the macro inserts " *
                  "the name automatically. Got $(length(ctor_args)) arguments. " *
                  "If you wrote the name explicitly as the second argument, " *
                  "remove it.")
        end

        name_str = String(key)
        new_call = Expr(:call, ctor, ctor_args[1], name_str, ctor_args[2], ctor_args[3])
        push!(rewritten_entries, Expr(:(=), key, new_call))
    end

    return esc(Expr(:tuple, rewritten_entries...))
end

end # module
