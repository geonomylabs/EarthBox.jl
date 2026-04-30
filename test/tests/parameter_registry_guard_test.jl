using Test
using EarthBox

import EarthBox.ParameterRegistry: get_eb_parameters

@testset "ParameterRegistry: NamedTuple key equals .name" begin
    for (k, v) in pairs(get_eb_parameters())
        @test String(k) == v.name
    end
end
