using Test
using Interferometers



@testset "Sliding modes: Minimum gain test" begin
    # close("all")
    τ = 1e-5; ϕ₀_amp = 1; f₀ = 3; f = 500; Δϕ_amp = 1 
    t  = 0:τ:0.005
    n = length(t)
    ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
    Δϕ = @. Δϕ_amp*sin(2π*f*t)
    ϕ = Δϕ .+ ϕ₀
    A1 = 1;    V1 = 1;
    arr_cos = A1*V1*cos.(ϕ); arr_sin = A1*V1*sin.(ϕ)

    gain_min = 2π*f₀*ϕ₀_amp + 2π*f*Δϕ_amp
    arr_gain = 0:100:30*gain_min
    e = range(start=0,stop=1e-1,length=11)
    max_error = zeros(size(arr_gain))
    _,ax1 = subplots()
    _,ax2 = subplots()
    for (j, e) in enumerate(e)
        for (i, gain) in enumerate(arr_gain)
            phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=π/2, solver=BS3)
            ϕc = -phase.phase
            max_error[i] = maximum(abs.(ϕ[Int(floor(n/2)):end] .+ ϕc[Int(floor(n/2)):end]))
            if e == 0.0 && gain == 3200 # 12600
                ax1.plot(t, -ϕc, label=L"$-\phi_c + \pi/2$")
                ax1.plot(t, ϕ, label="input")
                    ax1.set_ylabel("Phases")
                    ax1.set_title("gain=$(round(gain,digits=0)), dt=$(round(τ*1e6,digits=0)) us, e=$(round(e,digits=3))")
                    ax1.legend()
                ax11 = ax1.twinx()
                ax11.plot(t, ϕ.+ϕc, color="gray",label="error")
                    # ax11.set_ylim(-0.6,0.6); #ax11.set_xlim(0.003, 0.005)
                    ax11.legend()
            end

        end
        ax2.plot(arr_gain, max_error, label="$(round(e,digits=3))")
            ax2.set_title("dt=$(round(τ*1e6,digits=0)) us")

    end
    ax2.vlines(gain_min,0,2,linewidth=1,color="red",label="minimum gain")
        ax2.set_ylabel("Absolute error [rad]"); ax2.set_xlabel("Gain");
        ax2.set_yscale("log"); ax2.set_xscale("log")
        ax2.legend()
        ax2.grid(which="minor")
        ax2.set_ylim(0,3)
    
end





@testset "Sliding modes: sigmoid error and gain RK4" begin
    τ = 1e-6; ϕ₀_amp = 1; f₀ = 3; f = 500; Δϕ_amp = 1 
    t  = 0:τ:0.005
    n = length(t)
    ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
    Δϕ = @. Δϕ_amp*sin(2π*f*t)
    ϕ = Δϕ .+ ϕ₀
    A1 = 1;    V1 = 1;
    arr_cos = A1*V1*cos.(ϕ); arr_sin = A1*V1*sin.(ϕ)

    gain_min = 2π*f₀*ϕ₀_amp + 2π*f*Δϕ_amp

    arr_gain = 0:200:30*gain_min
    e = range(start=0,stop=1e-1,length=10)

    max_error = zeros(size(arr_gain))
    min_error_e = zeros(size(e))
    min_gain_e = zeros(size(e))
    for (j, e) in enumerate(e)
        for (i, gain) in enumerate(arr_gain)

            phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=0, solver=RK4)
            
            ϕc = -phase.phase
            max_error[i] = maximum(abs.(ϕ[Int(floor(n/2)):end] .+ ϕc[Int(floor(n/2)):end]))
            if e == -0.1 && gain == -321800
                figure(1)
                plot(t, -ϕc, label=L"$-\phi_c + \pi/2$")
                plot(t, ϕ, label="input")
                    ylabel("Phases"); title("Gain=$(round(gain,digits=0)), dt=$(round(τ*1e6,digits=0)) us"); legend()
                ax = twinx()
                plot(t, ϕ.+ϕc, color="gray",label="error")
                    legend(); ylim(-0.6,0.6); xlim(0.003, 0.005)
            end

        end
        figure(2)
        plot(arr_gain, max_error, label="$(round(e,digits=3))")
            title("dt=$(round(τ*1e6,digits=0)) us, sigm = $(e)")
        
        index = argmin(max_error)
        min_error_e[j] = max_error[index]
        min_gain_e[j] = arr_gain[index]

    end
    figure(2)
    vlines(gain_min,0,2,linewidth=1,color="red",label="minimum gain")
        ylabel("Absolute error [rad]"); xlabel("Gain"); yscale("log"); xscale("log")
        legend()
        grid(which="minor")

    figure()
    plot(e, min_error_e)
        xlabel(L"sigmoid factor $\varepsilon$")
        ylabel("Minimum absolute error [rad]")#, color=custom_plot_colors[1])
    ax = twinx()
    ax.plot(e, min_gain_e, color="tab:red")#, color=custom_plot_colors[2])
        ylabel("Gain for minimum error", color="tab:red")#, color=custom_plot_colors[2])
        # gca().spines["right"].set_visible(true)
    
end








@testset "Sliding modes: Thesis Felao 1 " begin

    """ Reproducing Fig 16 """
    function simulate_state(τ, tf, ϕ₀_amp, f₀, Δϕ_amp, f, gain, phase_ic)
        A1 = 1
        V1 = 1
        t  = 0:τ:tf
        ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t+phase_ic)
        Δϕ = @. Δϕ_amp*sin(2π*f*t)
        ϕ = Δϕ .+ ϕ₀
        arr_cos = A1*V1*cos.(ϕ)
        arr_sin = A1*V1*sin.(ϕ)

        gain_min = 2π*f₀*ϕ₀_amp + 2π*f*Δϕ_amp

        phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=0, ic=π/2, solver=BS3)
        ϕc = -phase.phase .+ π/2        # ϕ + ϕ₀ + ϕc = π/2

        @test(maximum(abs.(ϕ .+ ϕc)/π) ≈ 1/2, atol=2e-1)
        #TODO: check x1 and f convergence to x1=0 and f = -1 or +1
        _,ax = subplots(4,1,figsize=(5,5))
            ax[1].plot(t, -ϕc, label="-control")
            ax[1].plot(t, ϕ, label="input")
                ax[1].legend(loc="right"); ax[1].set_ylabel("Phases")
            ax[2].plot(t, -ϕc .- ϕ .+ π/2)
                ax[2].set_ylabel("error")    
            ax[3].plot(t, A1*V1*cos.(ϕ), label="open")
            ax[3].plot(t, A1*V1*cos.(ϕ .+ ϕc), label="closed")
                ax[3].legend(loc="right");ax[3].set_ylabel("x1"); ax[3].set_ylim(-1.1,1.1)
            ax[4].plot(t, -A1*V1*sin.(ϕ), label="open")
            ax[4].plot(t, -A1*V1*sin.(ϕ .+ ϕc), label="closed")
                ax[4].legend(loc="right"); ax[4].set_ylabel("f"); ax[4].set_ylim(-1.1,1.1)

    end
    
    tf = 1; gain = 3200
    ϕ₀_amp = 1; f₀ = 3
    Δϕ_amp = 1; f  = 500
    phase_ic = 0
    for τ in [1e-4 1e-5 1e-6]
        simulate_state(τ, tf, ϕ₀_amp, f₀, Δϕ_amp, f, gain, phase_ic)
        suptitle("period= $(τ*1e6) us")
    end


end





@testset "Sliding modes: Fig 18" begin
    
    τ = 1e-6
    t  = 0:τ:0.01
    ϕ₀_amp = 0; f₀ = 3; Δϕ_amp = 1; f = 500
    ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
    Δϕ = @. Δϕ_amp*sin(2π*f*t)
    ϕ = Δϕ .+ ϕ₀
    A1 = 1;    V1 = 1;  arr_cos = A1*V1*cos.(ϕ); arr_sin = A1*V1*sin.(ϕ)

    gain_min = 2π*f₀*ϕ₀_amp + 2π*f*Δϕ_amp
    gain = 3200;    ic = 0.7π
    
    phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=0, ic=ic, solver=BS3)
    ϕc = -phase.phase

    figure()
    plot(t, -ϕc, label=L"$-\phi_c + \pi/2$")
    plot(t, ϕ, label="input")
        legend()
        ylabel("Phases")
        title("Initial condition = \${$(ic/π)} \\pi\$")
    
    p1 = 500
    max_error = maximum(abs.(ϕ[p1:end] .+ ϕc[p1:end]))
    @test(max_error ≈ 0.0, atol=1e-1)

end






@testset "Sliding modes: offset due to SMC initial condition" begin
    
    ic = -2π:π/10:2π
    offset = zeros(size(ic))
    for (i,ic) in enumerate(ic)
        τ = 1e-6
        t  = 0:τ:0.01
        f₀ = 3
        f = 500
        ϕ₀_amp = 0
        Δϕ_amp = 1
        ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
        Δϕ = @. Δϕ_amp*sin(2π*f*t)
        ϕ = Δϕ .+ ϕ₀
        A1 = 1;    V1 = 1
        arr_cos = A1*V1*cos.(ϕ)
        arr_sin = A1*V1*sin.(ϕ)

        gain_min = 2π*f₀*ϕ₀_amp + 2π*f*Δϕ_amp
        gain = 3200
        phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=0, ic=ic)
        ϕc = -phase.phase

        offset[i] = phase.offset
        # TODO: implement test to check equilibrium points

    end
    figure()
        plot(ic/π, offset/π, ".")
        ylabel(L"offset [$\pi$ rad]")
        xlabel(L"Initial control condition [$\pi$ rad]")

end






@testset "Slinding modes boundary layer, thesis Felao 1" begin

    τ  = 10e-6
    t  = 0:τ:1/3
    ϕ₀_amp = 1
    Δϕ_amp = 1
    A1 = 1
    V1 = 1
    ϕ₀ = @. ϕ₀_amp*sin(2π*3*t)
    Δϕ = @. Δϕ_amp*sin(2π*500*t)

    C  = [3160; 3180; 3200; 4000:1000:15000]

    # getting the control signal by simulation to calculate the chattering
    arr_cos = cos.(Δϕ.+ϕ₀)
    arr_sin = sin.(Δϕ.+ϕ₀)
    ϕc = phase_highgain(arr_cos,arr_sin, τ, C[1], e=1e-6)
    ϕc = ϕc.phase

    # TODO: Chattering equation does not provide the thesis result
    chattering = maximum(-ϕc.+Δϕ.+ϕ₀)-minimum(-ϕc.+Δϕ.+ϕ₀)
    round.(chattering, digits=4)
    
    # figure()
    #     plot(t, ϕc,alpha=0.5)
    #     plot(t, Δϕ.+ϕ₀)
    # figure()
    #     plot(t, +ϕc.+Δϕ.+ϕ₀ .-π/2)


    # TODO: boundary layer equation does not provide the thesis result
    boundary_layer = @. 2*A1*V1*τ*(ϕ₀_amp + Δϕ_amp + C)
    round.(boundary_layer, digits=4)

    # figure()
    #     plot(t,arr_cos)
    # figure()
    #     plot(t, Δϕ)
    #     plot(t, ϕ₀)

end

