using Test
using EarthBox

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays

@testset "ArrayRegistry: NamedTuple key equals .name" begin
    for (k, v) in pairs(get_eb_arrays())
        @test String(k) == v.name
    end
end
