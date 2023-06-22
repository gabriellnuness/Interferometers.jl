"""
Simulation of a Michelson interferometer
"""
function sagnac(t, τ, ϕs, ϕm)
    0.5 + 0.5*cos(ϕs(t) + ϕm(t) - ϕm(t-τ))
end