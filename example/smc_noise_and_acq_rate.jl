# Simulate the figure from Felao thesis in which different acquisition rates
# are simulated, showing a smaller chattering to higher rates.
# But, what is there is white noise increasing together with the acq rate?

using Interferometers
using Statistics
using FFTW
using PyPlot

function white_noise(Δf=1)
    kb = 1.380649*10^(-23)
    T = 300 # ≈25°C
    R = 50
    voltage_rms = sqrt(4*kb*R*T*Δf)
    return voltage_rms
end

f = range(1, 1e6, 1000)
white_noise.(f)

figure()
plot(f, white_noise.(f))
    xlabel("Frequency")
    ylabel("RMS voltage noise")


# Firstly, create a random noise with the rms value
# the rms value of a signal with zero mean is the standard deviation
function make_signal()
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
close("all")
make_signal()




# testing different acquisition rate to increasing white noise
ϕ₀_amp = 1; f₀ = 3; f = 500; Δϕ_amp = 1 
τ = 2e-6
t  = 0:τ:0.005
ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
Δϕ = @. Δϕ_amp*sin(2π*f*t)
ϕ = Δϕ .+ ϕ₀
A1 = 1;    V1 = 1;
arr_cos = A1*V1*cos.(ϕ); arr_sin = A1*V1*sin.(ϕ)
gain = 2π*f₀*ϕ₀_amp + 2π*f*Δϕ_amp
e = 0

phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, ic=π/2, solver=BS3)
ϕc = -phase.phase

figure()
    plot(t, -ϕc, label=L"$-\phi_c + \pi/2$")
    plot(t, ϕ, label="input")
        ylabel("Phases")
        title("gain=$(round(gain,digits=0)), dt=$(round(τ*1e6,digits=0)) us, e=$(round(e,digits=3))")
        legend()
    ax1 = twinx()
    ax1.plot(t, ϕ.+ϕc, color="gray", label="error")
    ax1.legend()
