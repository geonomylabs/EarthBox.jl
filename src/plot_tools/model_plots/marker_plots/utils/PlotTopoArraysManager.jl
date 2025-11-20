module PlotTopoArraysManager

mutable struct PlotTopoArrays
    topox::Vector{Float64}
    topoy::Vector{Float64}
    
    function PlotTopoArrays()
        topox = zeros(Float64, 1)
        topoy = zeros(Float64, 1)
        new(topox, topoy)
    end
end

end # module
