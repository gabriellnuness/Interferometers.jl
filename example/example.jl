using PyPlot
using Interferometers

λ = 1310e-9
Δλ = 50e-9
ΔL = 100e-6

simout = fringes(λ, Δλ, ΔL, 1e-8)

plot(simout.length*1e6, simout.interferogram)
    xlabel(L"Optical path length ($\mu$m)")
    ylabel("Interference intensity")
    title("Low coherence interferogram")
