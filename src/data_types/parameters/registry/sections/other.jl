function get_other_parameters()::NamedTuple
    return @params (
        #*****************
        # Other Parameters
        #*****************
        iuse_mumps = ParameterInt(
            0, "None", "Flag to use MUMPS solver: 0 off; 1 on"),
    )
end
