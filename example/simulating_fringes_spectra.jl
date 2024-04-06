# Run this script in the main Julia Environment
# ]activate
using PyPlot
using Interferometers

# Simulating spatial interference fringes
λ = 1550e-9
Δλ = 150e-9
ΔL = 50e-6
simout = fringes(λ, Δλ, ΔL, Square)

figure()
plot(simout.length*1e6, simout.interferogram)
    xlabel(L"Optical path length ($\mu$m)")
    ylabel("Interference intensity")
    title("Low coherence interferogram")

simout = fringes(λ, Δλ, ΔL)

figure()
plot(simout.length*1e6, simout.interferogram)
    xlabel(L"Optical path length ($\mu$m)")
    ylabel("Interference intensity")
    title("Low coherence interferogram")


""" simulating different bandwidths """
λ = 1550e-9
Δλ = 40e-9
ΔL = 5e-6
simout = fringes(λ, Δλ, ΔL)
Δλ = 300e-9
simout1 = fringes(λ, Δλ, ΔL)

figure()
plot(simout.length*1e6, simout.interferogram)
plot(simout1.length*1e6, simout1.interferogram)
    xlabel(L"Optical path length ($\mu$m)")
    ylabel("Interference intensity")
    title("Low coherence interferogram\nscale factor change due to bandiwdth")
    legend([L"$\Delta\lambda$ = 40 nm",L"$\Delta\lambda$ = 300 nm"])
