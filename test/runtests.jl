using Interferometers
using Test
using PyPlot
using Printf





@testset "Simulation" begin

    @test michelson(0, 1, 1) == 2
    @test michelson(0.0, 1.0, 1.0) == 2
    @test michelson(deg2rad(90), 1, 1) == 1
    @test michelson(deg2rad(90), 0, 0) == 0

end




@testset "Quadrature" begin

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

    @test(rad2deg(phase) ≈ ϕ1-ϕ2 - 90, atol=1e-10)
    @test(A1/A2 ≈ gain_ratio, atol=1e-10)
    @test(offset_1 ≈ A1, atol=1e-10)
    @test(offset_2 ≈ A2, atol=1e-10)





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
    
    @test(rad2deg(phase) ≈ ϕ2-ϕ1 - 90, atol=1e-10) # inverted signals
    @test(A2/A1 ≈ gain_ratio, atol=1e-10)
    @test(offset_1 ≈ A2, atol=1e-10)
    @test(offset_2 ≈ A1, atol=1e-10)

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






@testset "Arctangent phase retrieval" begin

    # Testing with pure sine and cosine with pure frequency
    f = 60
    τ = 1/f
    t = 0 : τ/1000 : τ*2

    freq_mod = 2*f
    amp_mod = 10π
    dc_mod = 0.0*π

    Δϕ = @. amp_mod*sin(2π*freq_mod*t) + dc_mod
    signal_cos = cos.(Δϕ)
    signal_sin = sin.(Δϕ)

    # arctangent test    
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # test modulation signal retrieval
    @test maximum((Δϕ - phase)/Δϕ) <= 1e-10
    @test dc_mod - phase_offset <= 1e-10


    # Plot arctangent method
    _, ax = subplots(3,1)
    ax[1].plot(t, signal_cos)
        ax[1].set_ylabel("signal_cos")
    ax[2].plot(t, signal_sin)
        ax[2].set_ylabel("signal_sin")
    ax[3].plot(t, Δϕ, linewidth=5,color="black",alpha=0.2)
    ax[3].plot(t, phase)
        ax[3].set_ylabel("Phase retrieved") 
        str = @sprintf("phase offset:%.3fπ", (phase_offset/π))
        ax[3].legend(["Δϕ",str])
    suptitle("arc tangent method")




    # dc offset
    dc_mod = 0.1*π

    Δϕ = @. amp_mod*sin(2π*freq_mod*t) + dc_mod
    signal_cos = cos.(Δϕ)
    signal_sin = sin.(Δϕ)

    # arctangent test    
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # test modulation signal retrieval
    @test maximum((Δϕ - phase)/Δϕ) <= 1e-10
    @test dc_mod - phase_offset <= 1e-10



end

@testset "Sliding modes phase retrieval" begin

    # Testing with pure sine and cosine with pure frequency
    f = 60
    τ = 1/f
    t = 0 : τ/1000 : τ*2

    freq_mod = 2*f
    amp_mod = 10π
    dc_mod = 0.0*π

    Δϕ = @. amp_mod*sin(2π*freq_mod*t) + dc_mod
    signal_cos = cos.(Δϕ)
    signal_sin = sin.(Δϕ)
    
    # High gain test    
    dt = t[2]-t[1]
    gain = 4e4
    sigmoid_factor = 1e-1 # 0 = signal function
    (phase, phase_offset) = phase_highgain(signal_cos, signal_sin, dt, gain, sigmoid_factor)

    # test modulation signal retrieval
    @test maximum((Δϕ - phase)/Δϕ) <= 100e-6 
    @test(dc_mod ≈ phase_offset, atol=1e-2)


    # Plot sliding modes method
    _, ax = subplots(3,1)
    ax[1].plot(t, signal_cos)
        ax[1].set_ylabel("signal_cos")
    ax[2].plot(t, signal_sin)
        ax[2].set_ylabel("signal_sin")
    ax[3].plot(t, Δϕ, linewidth=5,color="black",alpha=0.2)
    ax[3].plot(t, (phase .- phase_offset))
        ax[3].set_ylabel("Phase retrieved") 
        str = @sprintf("phase offset:%.3fπ", (phase_offset/π))
        ax[3].legend(["Δϕ",str])
    suptitle("sliding modes method")
    

    
    # Testing with pure sine and cosine with pure frequency modulation dc
    dc_mod = 0.1*π

    Δϕ = @. amp_mod*sin(2π*freq_mod*t) + dc_mod
    signal_cos = cos.(Δϕ)
    signal_sin = sin.(Δϕ)
    
    # High gain test    
    dt = t[2]-t[1]
    gain = 4e4
    sigmoid_factor = 1e-1 # 0 = signal function
    (phase, phase_offset) = phase_highgain(signal_cos, signal_sin, dt, gain, sigmoid_factor)

    # test modulation signal retrieval
    @test maximum((Δϕ - phase)/Δϕ) <= 100e-6 
    @test(dc_mod ≈ phase_offset, atol=1e-2)


end 
