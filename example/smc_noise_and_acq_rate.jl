# Simulate the figure from Felao thesis in which different acquisition rates
# are simulated, showing a smaller chattering to higher rates.
# But, what is there is white noise increasing together with the acq rate?

using Interferometers
using Statistics
using FFTW
using PyPlot
using DSP

"""
Johnson noise formula implementation
"""
function white_noise(Δf=1)
    kb = 1.380649*10^(-23)
    T = 300 # ≈25°C
    R = 50
    voltage_rms = sqrt(4*kb*R*T*Δf)
    return voltage_rms
end

f = range(1, 1e6, 10000)
white_noise.(f)

figure()
plot(f, white_noise.(f))
    xlabel("Frequency")
    ylabel("RMS voltage noise")
    xscale("log")



function make_signal()
    # Firstly, create a random noise with the rms value
    # the rms value of a signal with zero mean is the standard deviation

    fs = 1_000
    t_max = 100
    t_arr = range(0, t_max, t_max*fs-1)
    N = length(t_arr)

    # signal in time
    signal = white_noise(fs)*randn(length(t_arr))
    println("std=$(std(signal))\nmean=$(mean(signal))")

    # fft
    spec_amplitude = 2/N*abs.(fftshift(fft(signal)))
    spec_frequency = fftshift(fftfreq(N, fs))

    _,ax = subplots(2,1)
    ax[1].plot(signal)
    ax[2].plot(spec_frequency, spec_amplitude)
        ax[2].set_xlim(0, fs/2)
        # ax[2].set_yscale("log")
end
make_signal()



function make_noisy_signal()

    τ = 1e-6
    ϕ₀_amp = 1; f₀ = 3; f = 500; Δϕ_amp = 1 
    A1 = 1;    V1 = 1;
    t  = 0:τ:0.1
    ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
    Δϕ = @. Δϕ_amp*sin(2π*f*t)
    ϕ = Δϕ .+ ϕ₀
    arr_cos_noise = 5e3*white_noise(1/τ)*randn(length(ϕ))
    arr_sin_noise = 5e3*white_noise(1/τ)*randn(length(ϕ))
    arr_cos = A1*V1*cos.(ϕ)
    arr_sin = A1*V1*sin.(ϕ)
    
    figure()
    plot(arr_cos)
    plot(arr_cos + arr_cos_noise)

    gain = 3200
    e = 0

    phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=π/2, solver=BS3)
    ϕc = -phase.phase

    arr_cos = A1*V1*cos.(ϕ) .+ arr_cos_noise
    arr_sin = A1*V1*sin.(ϕ) .+ arr_sin_noise
    phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=π/2, solver=BS3)
    ϕc_noisy = -phase.phase
    
    figure()
    plot(t, A1*V1*cos.(ϕ.+ϕc_noisy.+π/2),
        linewidth=0.4, label="w/ noise")
    plot(t, A1*V1*cos.(ϕ.+ϕc.+π/2),
        linewidth=0.4, label="dt=$(round(τ*1e6,digits=0)) us")
    # ylabel("A₁V₁cos(Δϕ+ϕ₀+ϕc)")
    ylabel(L"$A_1 V_1 \cos(\Delta\phi+\phi_0+\phi_c)$")
    ylim(-1,1)
    grid()
    legend()

end
make_noisy_signal()


function make_quantized_signal()

    τ = 1e-6
    ϕ₀_amp = 1; f₀ = 3; f = 500; Δϕ_amp = 1 
    A1 = 1;    V1 = 1;
    t  = 0:τ:0.01
    ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
    Δϕ = @. Δϕ_amp*sin(2π*f*t)
    ϕ = Δϕ .+ ϕ₀

    arr_cos = A1*V1*cos.(ϕ)
    arr_cos_quant = round.(A1*V1*cos.(ϕ), digits=2)
    
    figure()
    plot(arr_cos, label="perfect")
    plot(arr_cos_quant, label="quantized")
        legend()
end
make_quantized_signal()



function make_noisy_quantized_signal()

    τ = 1e-6
    ϕ₀_amp = 1; f₀ = 3; f = 500; Δϕ_amp = 1 
    A1 = 1;    V1 = 1;
    t  = 0:τ:0.01
    ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
    Δϕ = @. Δϕ_amp*sin(2π*f*t)
    ϕ = Δϕ .+ ϕ₀

    arr_cos_noise = 5e3*white_noise(1/τ)*randn(length(ϕ))
    arr_sin_noise = 5e3*white_noise(1/τ)*randn(length(ϕ))
    arr_cos = A1*V1*cos.(ϕ)
    arr_sin = A1*V1*sin.(ϕ)
    arr_cos_quant = round.(A1*V1*cos.(ϕ), digits=2)
    arr_sin_quant = round.(A1*V1*sin.(ϕ), digits=2)
    fcut = 5000
    arr_cos_quant_filt = filtfilt(digitalfilter(Lowpass(fcut, fs=1/τ), Butterworth(2)), arr_cos_quant)
    arr_sin_quant_filt = filtfilt(digitalfilter(Lowpass(fcut, fs=1/τ), Butterworth(2)), arr_sin_quant)
    
    
    figure()
    plot(arr_cos + arr_cos_noise, label="noisy")
    plot(arr_cos_quant, label="quantized")
    plot(arr_cos_quant_filt, label="quantized filtered")
    plot(arr_cos, label="perfect")
        legend()

    gain = 3200
    e = 0

    phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=π/2, solver=BS3)
    ϕc = -phase.phase

    arr_cos = A1*V1*cos.(ϕ) .+ arr_cos_noise
    arr_sin = A1*V1*sin.(ϕ) .+ arr_sin_noise
    phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=π/2, solver=BS3)
    ϕc_noisy = -phase.phase
    
    phase = phase_highgain(arr_cos_quant, arr_sin_quant, τ, gain, e=e, ic=π/2, solver=BS3)
    ϕc_quant = -phase.phase

    phase = phase_highgain(arr_cos_quant_filt, arr_sin_quant_filt, τ, gain, e=e, ic=π/2, solver=BS3)
    ϕc_quant_filt = -phase.phase

    
    figure()
    plot(t, A1*V1*cos.(ϕ.+ϕc_noisy.+π/2),
        linewidth=0.4, label="noisy")
    plot(t, A1*V1*cos.(ϕ.+ϕc_quant.+π/2),
        linewidth=0.4, label="quantized")
    plot(t, A1*V1*cos.(ϕ.+ϕc_quant_filt.+π/2),
        linewidth=0.4, label="quantized filtered")
    plot(t, A1*V1*cos.(ϕ.+ϕc.+π/2),
        linewidth=0.4, label="dt=$(round(τ*1e6,digits=0)) us")
    # ylabel("A₁V₁cos(Δϕ+ϕ₀+ϕc)")
    ylabel(L"$A_1 V_1 \cos(\Delta\phi+\phi_0+\phi_c)$")
    ylim(-1,1)
    grid()
    legend()

end
make_noisy_quantized_signal()












function noise_test()
    # using ifog 3x3 signal "amplitude"
    # testing different acquisition rate to increasing white noise

    τ_arr = range(1/100, 1/1_000_000, 10000)
    peak2peak = zeros(size(τ_arr))
    peak2peak_noise = zeros(size(τ_arr))

    # figure()
    for (i, τ) in enumerate(τ_arr)
        ϕ₀_amp = 1; f₀ = 1;
        Δϕ_amp = 10; f = 3;
        A1 = 1;    V1 = 1;
        t  = 0:τ:0.5
        ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
        Δϕ = @. Δϕ_amp*sin(2π*f*t)
        ϕ = Δϕ .+ ϕ₀
        arr_cos = A1*V1*cos.(ϕ); arr_sin = A1*V1*sin.(ϕ)
        arr_cos_noise = 5e3*white_noise(1/τ)*randn(length(ϕ))
        arr_sin_noise = 5e3*white_noise(1/τ)*randn(length(ϕ))

        gain = 2π*f₀*ϕ₀_amp + 2π*f*Δϕ_amp
        gain = gain*1.05
        e = 0

        phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=π/2, solver=BS3)
        ϕc = -phase.phase
        peak2peak[i] = maximum(A1*V1*cos.(ϕ.+ϕc.+π/2))-minimum(A1*V1*cos.(ϕ.+ϕc.+π/2))

        arr_cos = A1*V1*cos.(ϕ) .+ arr_cos_noise
        arr_sin = A1*V1*sin.(ϕ) .+ arr_sin_noise
        phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=π/2, solver=BS3)
        ϕc_noise = -phase.phase
        peak2peak_noise[i] = maximum(A1*V1*cos.(ϕ.+ϕc_noise.+π/2))-minimum(A1*V1*cos.(ϕ.+ϕc_noise.+π/2))

    end

    f = 1 ./τ_arr
    figure()
    plot(f, peak2peak, label="clean")
    plot(f, peak2peak_noise, label="noisy")
        xlabel("Frequency [Hz]")
        ylabel("peak to peak [V]")
        xscale("log")

    ax = twinx()
    ax.plot(f, white_noise.(f), color="black")

end
noise_test()







function lissajous_test(phase_mod = 100)
    # checking how noisy can the signal be and still allow the identification by
    # the lissajous figure.

    t = 0:0.001:0.5

    # assigning input signals as cosine or sine 
    Δϕ_amp = deg2rad(phase_mod)
    f = 3
    A1 = 1;    V1 = 1;
    A2 = 2;    V2 = 0.8
    dc1 = 1; dc2 = 2
    Δϕ = @. Δϕ_amp*sin(2π*f*t)
    ϕ = Δϕ
    ϕ1 = deg2rad(-120)
    ϕ2 = deg2rad(120)
    signal_1 = dc1 .+ A1*V1*cos.(ϕ .+ ϕ1)
    signal_2 = dc2 .+ A2*V2*cos.(ϕ .+ ϕ2)
    
    # correcting representation of sine and cosine    
    (phase, gain_ratio, offset_1, offset_2) = quadrature_fit(signal_1, signal_2)

    (signal_cos, signal_sin) = quadrature_set(signal_1, signal_2, phase, gain_ratio, offset_1, offset_2)
    
    # Optional plots
    figure()
    plot(signal_1, signal_2) # original signal   
    plot(signal_cos, signal_sin) # original signal   
        axis("equal"); title("Lissajous")
        legend(["Original", "In quadrature"])
        ylim(-5,5)
        xlim(-5,5)
    figure()
        plot(signal_1)
        plot(signal_2)

    println("phase=$(rad2deg(phase))\nratio=$(gain_ratio)\ndc1=$offset_1)\ndc2=$(offset_2)")


end
lissajous_test()
lissajous_test(10)