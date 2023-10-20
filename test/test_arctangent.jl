using Test
using Interferometers
using DSP


@testset "atan with offset: Lemes 1" begin

    x = 4π; f = 1e3; θ = 0 ; ϕ₀ = 3π/4
    f_sample = 1249.5*1e3; t_final = 1e-3

    t = 0 : 1/f_sample : t_final
    ϕ = @. ϕ₀ + x*sin(2π*f*t + θ) 
    
    signal_cos = cos.(ϕ)
    signal_sin = sin.(ϕ)

    close("all")
    figure(); plot(t, signal_cos); plot(t, signal_sin)

    # plotting phase wih atan and atan2
    figure(); plot(t, atan.(signal_sin./signal_cos)./π)
        legend(["atan"]); ylabel(L"phase [$\pi$ rad]"); ylim(-1,1)
    figure(); plot(t, atan.(signal_sin, signal_cos)./π)
        legend(["atan2"]); ylabel(L"phase [$\pi$ rad]")
    
    # arctangent Lemes
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # unwrap function
    phase_raw = @. atan(signal_sin, signal_cos)
    unwrapped_phase = unwrap(phase_raw)

    figure()
    plot(t, ϕ./π, color="black", alpha=.2, linewidth=4)
    plot(t, unwrapped_phase./π)
    plot(t, (phase)./π, "--")
        legend(["input", "atan2", "atan"])
        ylabel(L"phase [$\pi$ rad]")

end

@testset "atan with offset: Lemes 2" begin
    # first sample of the simulated phase is > 2π

    x = 4π; f = 1e3; θ = deg2rad(80) ; ϕ₀ = 3π/4
    f_sample = 1249.5*1e3; t_final = 2e-3

    t = 0 : 1/f_sample : t_final
    ϕ = @. ϕ₀ + x*sin(2π*f*t + θ) 
    
    signal_cos = cos.(ϕ)
    signal_sin = sin.(ϕ)

    close("all")
    figure(); plot(t, signal_cos); plot(t, signal_sin)

    # plotting phase wih atan and atan2
    figure(); plot(t, atan.(signal_sin./signal_cos)./π)
        legend(["atan"]); ylabel(L"phase [$\pi$ rad]"); ylim(-1,1)
    figure(); plot(t, atan.(signal_sin, signal_cos)./π)
        legend(["atan2"]); ylabel(L"phase [$\pi$ rad]")
    
    # arctangent Lemes
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # unwrap function
    phase_raw = @. atan(signal_sin, signal_cos)
    unwrapped_phase = unwrap(phase_raw)

    med_atan = (maximum(phase) + minimum(phase))/2/π
    med_atan2 = (maximum(unwrapped_phase) + minimum(unwrapped_phase))/2/π

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


@testset "atan with offset: Lemes 3" begin
    # first sample of the simulated phase is > 2π

    x = 4π; f = 1e3; θ = deg2rad(80) ; ϕ₀ = 1.3π
    f_sample = 1249.5*1e3; t_final = 2e-3

    t = 0 : 1/f_sample : t_final
    ϕ = @. ϕ₀ + x*sin(2π*f*t + θ) 
    
    signal_cos = cos.(ϕ)
    signal_sin = sin.(ϕ)

    close("all")
    figure(); plot(t, signal_cos); plot(t, signal_sin)

    # plotting phase wih atan and atan2
    figure(); plot(t, atan.(signal_sin./signal_cos)./π)
        legend(["atan"]); ylabel(L"phase [$\pi$ rad]"); ylim(-1,1)
    figure(); plot(t, atan.(signal_sin, signal_cos)./π)
        legend(["atan2"]); ylabel(L"phase [$\pi$ rad]")
    
    # arctangent Lemes
    (phase, phase_offset) = phase_atan(signal_cos, signal_sin)

    # unwrap function
    phase_raw = @. atan(signal_sin, signal_cos)
    unwrapped_phase = unwrap(phase_raw)

    med_atan = (maximum(phase) + minimum(phase))/2/π
    med_atan2 = (maximum(unwrapped_phase) + minimum(unwrapped_phase))/2/π

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


@testset "atan static phase recover range: Lemes 4" begin


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

        med_atan = (maximum(phase) + minimum(phase))/2/π
        med_atan2 = (maximum(unwrapped_phase) + minimum(unwrapped_phase))/2/π

        ϕ_atan[i] = med_atan
        ϕ_atan2[i] = med_atan2

    end

    figure(); plot(ϕ₀./π, ϕ_atan,"."); xlabel(L"$\phi_{0}$"); ylabel(L"$\phi_{0}$ atan")
    figure(); plot(ϕ₀./π, ϕ_atan2,"."); xlabel(L"$\phi_{0}$"); ylabel(L"${\phi}_{0}$ atan2")


end