using Interferometers
using Test
using PyPlot
using Statistics

@testset "Simulation" begin

    @test michelson(0, 1, 1) == 2
    @test michelson(0.0, 1.0, 1.0) == 2
    @test michelson(deg2rad(90), 1, 1) == 1
    @test michelson(deg2rad(90), 0, 0) == 0

end

@testset "Quadrature" begin

    threshold = 1e-10
    t = 0:0.001:0.5
    dc1 = 1; A1 = 1;
    dc2 = 2; A2 = 2;
    
    # assigning input signals as cosine or sine 
    cosine = @. dc1 + A1*cos(2π*10*t + deg2rad(0))
    sine = @. dc2 + A2*sin(2π*10*t + deg2rad(0))
    δ = 10
    (signal_cos, signal_sin) = check_ellipse_rotation(cosine, sine, δ)
    @test cosine == signal_cos
    @test sine == signal_sin

    (signal_cos, signal_sin) = check_ellipse_rotation(sine, cosine, δ)
    @test cosine == signal_cos
    @test sine == signal_sin

    # fitting quadrature
    ϕ1 = 60
    ϕ2 = -60
    signal_1_orig = @. dc1 + A1*cos(2π*10*t + deg2rad(ϕ1))
    signal_2_orig = @. dc2 + A2*cos(2π*10*t + deg2rad(ϕ2))
    # correcting representation of sine and cosine
    (signal_1, signal_2) = check_ellipse_rotation(signal_1_orig, signal_2_orig, δ)
    
    (phase, gain_ratio, offset_1, offset_2) = quadrature_fit(signal_1, signal_2)

    (signal_cos, signal_sin) = quadrature_set(signal_1, signal_2, phase, gain_ratio, offset_1, offset_2)
       
    # Optional plots
    figure()
    plot(signal_1_orig, signal_2_orig) # original signal
    x = signal_1_orig[1]; x1 = signal_1_orig[10]
    y = signal_2_orig[1]; y1 = signal_2_orig[10]
    arrow(x, y, x1-x, y1-y, head_width=0.1,color="black",zorder=10,label="_nolegend_")
    plot(signal_cos, signal_sin) # fitted signal
    x = signal_cos[1]; x1 = signal_cos[10]
    y = signal_sin[1]; y1 = signal_sin[10]
    arrow(x, y, x1-x, y1-y, head_width=0.1,color="black",zorder=10,label="_nolegend_")
        axis("equal"); title("Lissajous")
        legend(["Original", "In quadrature"])
        title("Not inverted")

    @test rad2deg(phase) - (ϕ1-ϕ2 - 90) <= threshold
    @test abs(A1/A2 - gain_ratio) <= threshold
    @test abs(offset_1 - A1) <= threshold
    @test abs(offset_2 - A2) <= threshold

    # inverting signals
    ϕ1 = -60
    ϕ2 = 60
    signal_1_orig = @. dc1 + A1*cos(2π*10*t + deg2rad(ϕ1))
    signal_2_orig = @. dc2 + A2*cos(2π*10*t + deg2rad(ϕ2))
    δ = 10
    (signal_1, signal_2) = check_ellipse_rotation(signal_1_orig, signal_2_orig, δ)
    # fitting quadrature
    (phase, gain_ratio, offset_1, offset_2) = quadrature_fit(signal_1, signal_2)
    # setting quadrature
    (signal_cos, signal_sin) = quadrature_set(signal_1, signal_2, phase, gain_ratio, offset_1, offset_2)
    
    @test rad2deg(phase) - (ϕ2-ϕ1 - 90) <= threshold # inverted signals
    @test abs(A2/A1 - gain_ratio) <= threshold
    @test abs(offset_1 - A2) <= threshold
    @test abs(offset_2 - A1) <= threshold

    figure()
    plot(signal_1_orig, signal_2_orig) # original signal
    x = signal_1_orig[1]; x1 = signal_1_orig[10]
    y = signal_2_orig[1]; y1 = signal_2_orig[10]
    arrow(x, y, x1-x, y1-y, head_width=0.1,color="black",zorder=10,label="_nolegend_")
    plot(signal_cos, signal_sin) # fitted signal
    x = signal_cos[1]; x1 = signal_cos[10]
    y = signal_sin[1]; y1 = signal_sin[10]
    arrow(x, y, x1-x, y1-y, head_width=0.1,color="black",zorder=10,label="_nolegend_")
        axis("equal"); title("Lissajous")
        legend(["Original", "In quadrature"])
        title("Inverted")

end

        
