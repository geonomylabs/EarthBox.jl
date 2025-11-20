module CalculateMeanResiduals

import ..LevelManager: LevelData, LevelData2d

function calculate_scaled_and_mean_residuals!(
    ivcycle::Int64,
    stokesscale::Float64,
    continscale::Float64,
    level1_data::LevelData,
    ΔRxL::Array{Float64,3},
    ΔRyL::Array{Float64,3},
    ΔRzL::Array{Float64,3},
    ΔRcL::Array{Float64,3},
    resx::Vector{Float64},
    resy::Vector{Float64},
    resz::Vector{Float64},
    resc::Vector{Float64}
)::Tuple{Array{Float64,3}, Array{Float64,3}, Array{Float64,3}, Array{Float64,3}}
    ynum = level1_data.grid.parameters.geometry.ynum.value
    xnum = level1_data.grid.parameters.geometry.xnum.value
    znum = level1_data.grid.parameters.geometry.znum.value

    resx[ivcycle] = 0.0
    resy[ivcycle] = 0.0
    resz[ivcycle] = 0.0
    resc[ivcycle] = 0.0
    
    ΔRxL_scaled = ΔRxL ./ stokesscale
    ΔRyL_scaled = ΔRyL ./ stokesscale
    ΔRzL_scaled = ΔRzL ./ stokesscale
    ΔRcL_scaled = ΔRcL ./ continscale

    for k = 1:znum
        for j = 1:xnum
            for i = 1:ynum
                # x-Stokes
                if i > 1 && j > 1 && k > 1 && j < xnum
                    resx[ivcycle] += ΔRxL_scaled[i,j,k]^2
                end
                # y-Stokes
                if i > 1 && j > 1 && k > 1 && i < ynum
                    resy[ivcycle] += ΔRyL_scaled[i,j,k]^2
                end
                # z-Stokes
                if i > 1 && j > 1 && k > 1 && k < znum
                    resz[ivcycle] += ΔRzL_scaled[i,j,k]^2
                end
                # Continuity
                if j < xnum && i < ynum && k < znum
                    resc[ivcycle] += ΔRcL_scaled[i,j,k]^2
                end
            end
        end
    end
    resx[ivcycle] = log10(sqrt(resx[ivcycle] / ((ynum-1) * (xnum-2) * (znum-1))))
    resy[ivcycle] = log10(sqrt(resy[ivcycle] / ((ynum-2) * (xnum-1) * (znum-1))))
    resz[ivcycle] = log10(sqrt(resz[ivcycle] / ((ynum-1) * (xnum-1) * (znum-2))))
    resc[ivcycle] = log10(sqrt(resc[ivcycle] / ((ynum-1) * (xnum-1) * (znum-1))))
    return ΔRxL_scaled, ΔRyL_scaled, ΔRzL_scaled, ΔRcL_scaled
end

function calculate_scaled_and_mean_residuals!(
    ivcycle::Int64,
    stokesscale::Float64,
    continscale::Float64,
    level1_data::LevelData2d,
    ΔRxL::Array{Float64,2},
    ΔRyL::Array{Float64,2},
    ΔRcL::Array{Float64,2},
    resx::Vector{Float64},
    resy::Vector{Float64},
    resc::Vector{Float64}
)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}
    ynum = level1_data.grid.parameters.geometry.ynum.value
    xnum = level1_data.grid.parameters.geometry.xnum.value

    resx[ivcycle] = 0.0
    resy[ivcycle] = 0.0
    resc[ivcycle] = 0.0
    
    ΔRxL_scaled = ΔRxL ./ stokesscale
    ΔRyL_scaled = ΔRyL ./ stokesscale
    ΔRcL_scaled = ΔRcL ./ continscale

    for j = 1:xnum
        for i = 1:ynum
            # x-Stokes
            if i > 1 && j > 1 && j < xnum
                resx[ivcycle] += ΔRxL_scaled[i,j]^2
            end
            # y-Stokes
            if i > 1 && j > 1 && i < ynum
                resy[ivcycle] += ΔRyL_scaled[i,j]^2
            end
            # Continuity
            if j < xnum && i < ynum
                resc[ivcycle] += ΔRcL_scaled[i,j]^2
            end
        end
    end
    resx[ivcycle] = log10(sqrt(resx[ivcycle] / ((ynum-1) * (xnum-2))))
    resy[ivcycle] = log10(sqrt(resy[ivcycle] / ((ynum-2) * (xnum-1))))
    resc[ivcycle] = log10(sqrt(resc[ivcycle] / ((ynum-1) * (xnum-1))))
    return ΔRxL_scaled, ΔRyL_scaled, ΔRcL_scaled
end

end