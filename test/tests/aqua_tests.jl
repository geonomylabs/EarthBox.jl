using Aqua
using EarthBox

@testset "Aqua quality checks" begin
    Aqua.test_all(EarthBox; ambiguities = (recursive = false,))
end
