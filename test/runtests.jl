using Interferometers
using Test

@testset "Interferometers.jl" begin

    @test michelson(0, 1, 1) == 2
    @test michelson(0.0, 1.0, 1.0) == 2
    @test michelson(deg2rad(90), 1, 1) == 1
    @test michelson(deg2rad(90), 0, 0) == 0

end
