module SysTools

function get_username()::String
    return get(ENV, "USER", get(ENV, "USERNAME", "unknown"))
end

end # module