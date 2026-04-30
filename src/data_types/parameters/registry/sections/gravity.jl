function get_gravity_parameters()::NamedTuple
    return @params (
        gravity_x = ParameterFloat(
            0.0, "m/s/s", "x-component of gravity acceleration in m/s/s"),
        gravity_y = ParameterFloat(
            9.8, "m/s/s", "y-component of gravity acceleration in m/s/s"),
        turn_off_gravity_y = ParameterInt(
            0, "None", "Flag to turn off y-component of gravity: 0 = off; 1 = on"),
        nsteps_turn_off_gravity = ParameterInt(
            100, "None", "Number of steps to gradually turn off gravity"),
    )
end
