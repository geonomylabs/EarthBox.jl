module PressureUpdate

import ...LevelManager: LevelData, LevelData2d

function update_pressure!(
    prnorm::Float64,
    level_data::LevelData
)::Nothing
    pr = level_data.pr.array
    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value
    znum = level_data.grid.parameters.geometry.znum.value
    dp = prnorm - pr[1,1,1]
    for k = 1:znum-1
        for j = 1:xnum-1
            for i = 1:ynum-1
                pr[i,j,k] += dp
            end
        end
    end
    return nothing
end

function update_pressure!(
    prnorm::Float64,
    level_data::LevelData2d
)::Nothing
    pr = level_data.pr.array
    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value
    dp = prnorm - pr[1,1,1]
    for j = 1:xnum-1
        for i = 1:ynum-1
            pr[i,j] += dp
        end
    end
    return nothing
end

end # module