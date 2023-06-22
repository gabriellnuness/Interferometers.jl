using Interferometers
using Test
using PyPlot

@testset "Simulation" begin

    @test michelson(0, 1, 1) == 2
    @test michelson(0.0, 1.0, 1.0) == 2
    @test michelson(deg2rad(90), 1, 1) == 1
    @test michelson(deg2rad(90), 0, 0) == 0

end

@testset "Interferometric phase" begin

    t = 0:0.001:0.5
    ϕ1 = -60; dc1 = 1; A1 = 1;
    ϕ2 = 60;  dc2 = 2; A2 = 2;
    signal_1 = @. dc1 + A1*sin(2π*10*t + deg2rad(ϕ1))
    signal_2 = @. dc2 + A2*sin(2π*10*t + deg2rad(ϕ2))

    # fitting quadrature
    (phase, gain_ratio, offset_1, offset_2) = quadrature_fit(signal_1, signal_2)
    @test round(rad2deg(phase)) == ϕ2-ϕ1 - 90
    @test A1/A2 - gain_ratio <= 1e-4
    @test offset_1 - A1 <= 1e-4
    @test offset_2 - A2 <= 1e-4

    (signal_1_cos, signal_2_sin) = quadrature_set(signal_1, signal_2,phase, gain_ratio, offset_1, offset_2)

    # optional plots
    plot(signal_1,signal_2)
    plot(signal_1_cos,signal_2_sin)
        axis("equal"); title("Lissajous")
        legend(["Original", "In quadrature"])
end

        