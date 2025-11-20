module ResidualStructs

struct StokesResidualsInput
    xnum::Int
    ynum::Int
    xstp_b::Vector{Float64}
    ystp_b::Vector{Float64}
    xstp_vy::Vector{Float64}
    ystp_vx::Vector{Float64}
    bintern_zone::Vector{Int64}
    etas0::Matrix{Float64}
    etan0::Matrix{Float64}
    RX1::Matrix{Float64}
    RY1::Matrix{Float64}
    RC1::Matrix{Float64}
    vx1::Matrix{Float64}
    vy1::Matrix{Float64}
    pr1::Matrix{Float64}
end

end # module 