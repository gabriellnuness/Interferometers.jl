using Test
using Interferometers
using DSP

"""
Tests for arctangent phase unwrap method,
including exact values from Lemes' thesis.
"""
plot_flag = false

@testset "Arctangent Lemes 1: ϕ₀<π, θ=0" begin
    x = 4π; f = 1e3; θ = 0 ; ϕ₀ = 3π/4
    f_sample = 1249.5*1e3; t_final = 1e-3

    t = 0 : 1/f_sample : t_final
    ϕ = @. ϕ₀ + x*sin(2π*f*t + θ) 
    signal_cos = cos.(ϕ)
    signal_sin = sin.(ϕ)

    # arctangent Lemes
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # unwrap function
    phase_raw = @. atan(signal_sin, signal_cos)
    unwrapped_phase = unwrap(phase_raw)

    if plot_flag
        figure(); plot(t, signal_cos); plot(t, signal_sin)
        figure(); plot(t, atan.(signal_sin./signal_cos)./π)
            legend(["atan"]); ylabel(L"phase [$\pi$ rad]"); ylim(-1,1)
        figure(); plot(t, atan.(signal_sin, signal_cos)./π)
            legend(["atan2"]); ylabel(L"phase [$\pi$ rad]")
        figure()
            plot(t, ϕ./π, color="black", alpha=.2, linewidth=4)
            plot(t, unwrapped_phase./π)
            plot(t, (phase)./π, "--")
                legend(["input", "atan2", "atan"])
                ylabel(L"phase [$\pi$ rad]")
    end
end

@testset "Arctangent Lemes 2: ϕ₀<π, θ=80°" begin
    # first sample of the simulated phase is > 2π
    x = 4π; f = 1e3; θ = deg2rad(80) ; ϕ₀ = 3π/4
    f_sample = 1249.5*1e3; t_final = 2e-3

    t = 0 : 1/f_sample : t_final
    ϕ = @. ϕ₀ + x*sin(2π*f*t + θ) 
    signal_cos = cos.(ϕ)
    signal_sin = sin.(ϕ)

    # arctangent Lemes
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # unwrap function
    phase_raw = @. atan(signal_sin, signal_cos)
    unwrapped_phase = unwrap(phase_raw)

    if plot_flag
        figure(); plot(t, signal_cos); plot(t, signal_sin)
        figure(); plot(t, atan.(signal_sin./signal_cos)./π)
            legend(["atan"]); ylabel(L"phase [$\pi$ rad]"); ylim(-1,1)
        figure(); plot(t, atan.(signal_sin, signal_cos)./π)
            legend(["atan2"]); ylabel(L"phase [$\pi$ rad]")
        figure()
            plot(t, ϕ./π, color="black", alpha=.2, linewidth=4)
            plot(t, unwrapped_phase./π)
            plot(t, (phase)./π, "--")
                legend(["input", "atan2", "atan"])
                ylabel(L"phase [$\pi$ rad]")
                title("Without phase removal from atan2")
        figure()
            plot(t, ϕ./π, color="black", alpha=.2, linewidth=4)
            plot(t, (unwrapped_phase .- atan_phase_offset(unwrapped_phase)[1])./π)
            plot(t, (phase)./π, "--")
                legend(["input", "atan2", "atan"])
                ylabel(L"phase [$\pi$ rad]")
                title("Removed phase offset from atan2")
    end
end


@testset "Arctangent Lemes 3: ϕ₀>π, θ=80°" begin
    # first sample of the simulated phase is > 2π
    x = 4π; f = 1e3; θ = deg2rad(80) ; ϕ₀ = 1.3π
    f_sample = 1249.5*1e3; t_final = 2e-3

    t = 0 : 1/f_sample : t_final
    ϕ = @. ϕ₀ + x*sin(2π*f*t + θ) 
    signal_cos = cos.(ϕ)
    signal_sin = sin.(ϕ)

    # arctangent Lemes
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # unwrap function
    phase_raw = @. atan(signal_sin, signal_cos)
    unwrapped_phase = unwrap(phase_raw)
    unwrapped_phase = unwrapped_phase .- atan_phase_offset(unwrapped_phase)[1]

    if plot_flag
        figure(); plot(t, signal_cos); plot(t, signal_sin)
        figure(); plot(t, atan.(signal_sin./signal_cos)./π)
            legend(["atan"]); ylabel(L"phase [$\pi$ rad]"); ylim(-1,1)
        figure(); plot(t, atan.(signal_sin, signal_cos)./π)
            legend(["atan2"]); ylabel(L"phase [$\pi$ rad]")
        figure()
            plot(t, ϕ./π, color="black", alpha=.2, linewidth=4)
            plot(t, unwrapped_phase./π)
            plot(t, (phase)./π, "--")
                legend(["input", "atan2", "atan"])
                ylabel(L"phase [$\pi$ rad]")
                title("Without phase removal from atan2")
        figure()
            plot(t, ϕ./π, color="black", alpha=.2, linewidth=4)
            plot(t, (unwrapped_phase .- atan_phase_offset(unwrapped_phase)[1])./π)
            plot(t, (phase)./π, "--")
                legend(["input", "atan2", "atan"])
                ylabel(L"phase [$\pi$ rad]")
                title("Removed phase offset from atan2")
    end
end


@testset "Arctangent Lemes 4: -π <ϕ₀<π, θ=80°" begin
    ϕ₀ = -2π : 2π/100 : 2π
    ϕ_atan = zeros(length(ϕ₀))
    ϕ_atan2 = zeros(length(ϕ₀))
        
    x = 4π; f = 1e3; θ = deg2rad(80)
    f_sample = 1249.5*1e3; t_final = 2e-3
    t = 0 : 1/f_sample : t_final

    for (i, ϕ₀) in enumerate(ϕ₀)

        ϕ = @. ϕ₀ + x*sin(2π*f*t + θ) 
        signal_cos = cos.(ϕ)
        signal_sin = sin.(ϕ)
        # arctangent Lemes
        (phase, phase_offset) = phase_atan(signal_cos, signal_sin)
        # unwrap function
        phase_raw = @. atan(signal_sin, signal_cos)
        unwrapped_phase = unwrap(phase_raw)
        unwrapped_phase = unwrapped_phase .- atan_phase_offset(unwrapped_phase)[1]

        ϕ_atan[i] = (maximum(phase) + minimum(phase))/2/π
        ϕ_atan2[i] = (maximum(unwrapped_phase) + minimum(unwrapped_phase))/2/π

    end

    if plot_flag
        figure(); plot(ϕ₀./π, ϕ_atan,".",markersize=2)
            xlabel(L"$\phi_{0}$"); ylabel(L"$\phi_{0}$ atan")
        figure(); plot(ϕ₀./π, ϕ_atan2,".",markersize=2)
            xlabel(L"$\phi_{0}$"); ylabel(L"${\phi}_{0}$ atan2")
    end
end



@testset "Arctangent fake wrap - aliasing" begin
    # Testing with pure sine and cosine with pure frequency
    f = 10 /(2π); ϕ = 20
    
    # calculation of limits to unwrap correctly fs > 2*f*ϕ
    fs = ceil(2*f*ϕ)
    ϕ_max = fs/(2f)
    
    t = 0 : 1/fs : 2*(1/f)
    Δϕ = @.  ϕ*sin(2π*f*t)
    signal_cos = cos.(Δϕ); signal_sin = sin.(Δϕ)

    # arctangent 
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)
    # comparing with DSP unwrap algorithm
    phase_raw = @. atan(signal_sin, signal_cos)
    unwrapped_phase = unwrap(phase_raw)
    unwrapped_phase = unwrapped_phase  .- atan_phase_offset(unwrapped_phase)[1]

    if plot_flag
        figure(); plot(t, signal_cos,".-",markersize=2); ylabel("input\nsignal")
        figure(); plot(t, phase_raw,".-",markersize=2); ylabel("atan2")
        figure()
            plot(t, Δϕ,".-",markersize=2,linewidth=5,color="black",alpha=0.2)
            plot(t, phase,".-",markersize=2)
            plot(t, unwrapped_phase,".-",markersize=2)
                legend(["input","unwrap home","unwrap DSP"]); ylabel("Phase [rad]")
    end
end








@testset "Arctangent fake wrap (Lemes test)" begin
    # first sample of the simulated phase is > 2π
    x = 4π; f = 1e3; θ = deg2rad(80) ; ϕ₀ = 3π/4
    t_final = 2e-3
    
    # limiting sample frequency: Φ < fs⋅π/(2πf)
    fs_min = 2*x*f # atan from Lemes needs twice as fast sample rate

    f_sample = fs_min
    t = 0 : 1/f_sample : t_final
    ϕ = @. ϕ₀ + x*sin(2π*f*t + θ) 
    signal_cos = cos.(ϕ)
    signal_sin = sin.(ϕ)

    # arctangent Lemes
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # unwrap function
    phase_raw = @. atan(signal_sin, signal_cos)
    unwrapped_phase = unwrap(phase_raw)

    if plot_flag
        figure()
            plot(t, ϕ./π, color="black", alpha=.2, linewidth=4)
            plot(t, (unwrapped_phase .- atan_phase_offset(unwrapped_phase)[1])./π)
            plot(t, (phase)./π, "--")
                legend(["input", "atan2", "atan"])
                ylabel(L"phase [$\pi$ rad]")
                title("fs = $(round(fs_min/1000,digits=3)) kHz\nInstead of previous fs = 1.2495 MHz")
    end
end