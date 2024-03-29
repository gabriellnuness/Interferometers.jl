using Interferometers
using Test
using Printf
using DSP
using PyPlot


@testset "Simulation" begin

    @test michelson(0, 1, 1) == 2
    @test michelson(0.0, 1.0, 1.0) == 2
    @test michelson(deg2rad(90), 1, 1) == 1
    @test michelson(deg2rad(90), 0, 0) == 0


    R = .7
    λ = 1550e-9
    L = range(start=0, stop=5e-6, length=10000)
    spec = fabry_perot.(R, L, λ)
    close("all")
    plot(L*1e6, spec)
        ylim(-0.2,1.2); xlabel("L [m]"); ylabel("Transmission")


end

@testset "Integrators" begin
    include("test_integrators.jl")
end


@testset "Quadrature" begin

    t = 0:0.001:0.5
    dc1 = 1; A1 = 1;
    dc2 = 2; A2 = 2;
    
    # assigning input signals as cosine or sine 
    cosine = @. dc1 + A1*cos(2π*10*t + deg2rad(0))
    sine = @. dc2 + A2*sin(2π*10*t + deg2rad(0))
    δ = 10
    (signal_cos, signal_sin) = make_cos_first(cosine, sine, δ)
    @test cosine == signal_cos
    @test sine == signal_sin

    (signal_cos, signal_sin) = make_cos_first(sine, cosine, δ)
    @test cosine == signal_cos
    @test sine == signal_sin

    # fitting quadrature
    ϕ1 = 60
    ϕ2 = -60
    signal_1_orig = @. dc1 + A1*cos(2π*10*t + deg2rad(ϕ1))
    signal_2_orig = @. dc2 + A2*cos(2π*10*t + deg2rad(ϕ2))
    
    # correcting representation of sine and cosine
    (signal_1, signal_2) = make_cos_first(signal_1_orig, signal_2_orig, δ)
    
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



    # test named tuple as input of next function
    quad_tuple = quadrature_fit(signal_1, signal_2)
    (signal_cos_tuple, signal_sin_tuple) = quadrature_set(signal_1, signal_2, quad_tuple)

    @test signal_cos == signal_cos_tuple
    @test signal_sin == signal_sin_tuple
       


    # non pure sine and cosine test
    ϕ1 = -60
    ϕ2 = 60
    signal_1_orig = @. dc1 + A1*cos(2π*10*t + deg2rad(ϕ1))
    signal_2_orig = @. dc2 + A2*cos(2π*10*t + deg2rad(ϕ2))
    δ = 10
    (signal_1, signal_2) = make_cos_first(signal_1_orig, signal_2_orig, δ)
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
        str = @sprintf("phase offset:%.2fpi", (phase_offset/π))
        ax[3].legend([L"\Delta\phi",str])
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



@testset "Artangent general tests" begin
    include("test_arctangent.jl")
end





@testset "Sliding modes phase retrieval" begin
    function plot_highgain()
        # Plot sliding modes method
        _, ax = subplots(3,1, figsize=(5,5))
        ax[1].plot(t, signal_cos)
            ax[1].set_ylabel("signal_cos")
        ax[2].plot(t, signal_sin)
            ax[2].set_ylabel("signal_sin")
        ax[3].plot(t, Δϕ, linewidth=3,color="black",alpha=0.2)
        ax[3].plot(t, (highgain.phase))
            ax[3].set_ylabel("Phase retrieved") 
            str = @sprintf("phase offset:%.2f pi", (highgain.offset/π))
            ax[3].legend([L"$\Delta\phi$",str])
        suptitle("sliding modes method")
        # test modulation signal retrieval
    end

    """ Testing with pure sine and cosine with pure frequency """
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
    gain_min = amp_mod*2π*freq_mod
    gain = gain_min + 50
    sigmoid_factor = 0 # 0 = signal function
    highgain = phase_highgain(signal_cos, signal_sin, dt, gain,
                                solver=BS3, e=sigmoid_factor, ic=π/2)
    plot_highgain()
    
    @test maximum((Δϕ - highgain.phase)/Δϕ) <= 100e-6 
    @test(dc_mod ≈ highgain.offset, atol=1e-2)





    """ Testing with pure sine and cosine with pure frequency modulation dc """
    dc_mod = 0.1*π

    Δϕ = @. amp_mod*sin(2π*freq_mod*t) + dc_mod
    signal_cos = cos.(Δϕ)
    signal_sin = sin.(Δϕ)
    
    # High gain test    
    dt = t[2]-t[1]
    gain = 4e4
    sigmoid_factor = 1e-1 # 0 = signal function
    highgain = phase_highgain(signal_cos, signal_sin, dt, gain, e=sigmoid_factor,solver=BS3,ic=π/2)
    plot_highgain()

    # test modulation signal retrieval
    @test maximum((Δϕ - highgain.phase)/Δϕ) <= 100e-6 
    @test(dc_mod ≈ highgain.offset, atol=1e-2)
   
   
   
   
   
   
    """ Testing modulation dc from -2π to 2π """
    dc_mod = -2π : 2π/100 : 2π
    offset = zeros(length(dc_mod))
    for (i, dc_mod) in enumerate(dc_mod)
        Δϕ = @. amp_mod*sin(2π*freq_mod*t) + dc_mod
        signal_cos = cos.(Δϕ)
        signal_sin = sin.(Δϕ)
        
        # High gain test    
        dt = t[2]-t[1]
        gain = 4e4
        sigmoid_factor = 1e-1 # 0 = signal function
        highgain = phase_highgain(signal_cos, signal_sin, dt, gain, e=sigmoid_factor,solver=BS3,ic=π/2)
        # plot_highgain()
        offset[i] = highgain.offset
        
        # test modulation signal retrieval
        if  -π/2 < dc_mod < π/2 
            @test maximum((Δϕ - highgain.phase)/Δϕ) <= 100e-6 
            @test(dc_mod ≈ highgain.offset, atol=1e-2)
        end
    end

    figure()
    plot(dc_mod./π, cos.(dc_mod .+ π/2), color="black",alpha=0.5)
    plot(dc_mod./π, offset./π, ".",markersize=2)
        # ylim(-.6,.6)
        ylabel(L"SMC offset [$\pi$ rad]"); xlabel(L"input phase offset [$\pi$ rad]")
        title(L"Phase offset SMC recovery\\$\Delta\phi = \Phi \sin(2\pi{\cdot}f{\cdot}t)+\phi_{dc}$")
end 


@testset "Sliding modes general tests" begin
    include("test_sliding_modes.jl")
end
