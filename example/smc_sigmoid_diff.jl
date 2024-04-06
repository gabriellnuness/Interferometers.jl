# Script to simulate the effect of the acquisition rate and sigmoid factor 
# on the optimal gain
using Test
using Interferometers


dt = [10;5;1;0.5]
gain_euler = [6.54;12.8;62.77;124.03]
gain_bs3 = [3.98;8.64;41.25;82.46]
gain_rk4 = [5.43;10.4;50.2;98.88]
dt2 = range(start=dt[1], stop=dt[end], length=100)


figure()
plot(dt, gain_euler, label="euler")

sol = Interferometers.fit_power(dt, gain_euler)
plot(dt2, sol[2][1].*dt2.^(sol[2][2]),label="\$ y = {$(round(sol[2][1],digits=2))} \\cdot x ^{$(round(sol[2][2],digits=2))}\$")

plot(dt, gain_bs3, label="bs3")
sol = Interferometers.fit_power(dt, gain_bs3)
plot(dt2, sol[2][1].*dt2.^(sol[2][2]),label="\$ y = {$(round(sol[2][1],digits=2))} \\cdot x ^{$(round(sol[2][2],digits=2))}\$")

plot(dt, gain_rk4, label="rk4")
sol = Interferometers.fit_power(dt, gain_rk4)
plot(dt2, sol[2][1].*dt2.^(sol[2][2]),label="\$ y = {$(round(sol[2][1],digits=2))} \\cdot x ^{$(round(sol[2][2],digits=2))}\$")

ylabel("Optimal gain / 0.1 sigmoid factor")
xlabel(L"dt [$\mu$s]")
legend()
println(sol[2])





# generate multiple points to increase plot resolution
sample_rate = range(start=0.5e-6, stop=10e-6, length=20)
coef = zeros(size(sample_rate))
for  (sample,τ) in enumerate(sample_rate)
    ϕ₀_amp = 1; f₀ = 3; f = 500; Δϕ_amp = 1 
    t  = 0:τ:0.005
    n = length(t)
    ϕ₀ = @. ϕ₀_amp*sin(2π*f₀*t)
    Δϕ = @. Δϕ_amp*sin(2π*f*t)
    ϕ = Δϕ .+ ϕ₀
    A1 = 1;    V1 = 1;
    arr_cos = A1*V1*cos.(ϕ); arr_sin = A1*V1*sin.(ϕ)

    solver= Euler
    gain_min = 2π*f₀*ϕ₀_amp + 2π*f*Δϕ_amp

    arr_gain = 0:200:100*gain_min
    e = range(start=0,stop=0.5e-1,length=11) # stop=1e-1

    max_error = zeros(size(arr_gain))
    min_error_e = zeros(size(e))
    min_gain_e = zeros(size(e))
    for (j, e) in enumerate(e)
        for (i, gain) in enumerate(arr_gain)

            phase = phase_highgain(arr_cos, arr_sin, τ, gain, e=e, solver=solver)
            
            ϕc = -phase.phase
            max_error[i] = maximum(abs.(ϕ[Int(floor(n/2)):end] .+ ϕc[Int(floor(n/2)):end]))

        end

        index = argmin(max_error)
        min_error_e[j] = max_error[index]
        min_gain_e[j] = arr_gain[index]

    end

    min_gain_e[1] = gain_min
    fit = Interferometers.fit_line(e, min_gain_e./min_gain_e[1])
    coef[sample] = fit[2][1]/10
    println(τ*1e6)
end 
figure()
plot(sample_rate.*1e6, coef)
    xlabel(L"dt [$\mu$s]")
    ylabel("Optimal gain / 0.1 sigmoid factor")
sol = Interferometers.fit_power(sample_rate.*1e6, coef)
plot(sample_rate.*1e6, sol[1])
legend(["data","\$ y = {$(round(sol[2][1],digits=2))} \\cdot x ^{$(round(sol[2][2],digits=2))}\$"])
println(sol[2])
