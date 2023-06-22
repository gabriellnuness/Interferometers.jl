"""
Simulation of a Michelson interferometer
"""
function michelson(ϕ,A,V)
    A + A*V*cos(ϕ)
end
