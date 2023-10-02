module Interferometers

using LinearAlgebra

# Simulation functions
export michelson
export sagnac
export fringes
# Phase recovery functions
export phase_atan
export phase_highgain
export make_cos_first
export quadrature_fit
export quadrature_set


include("simulation_models.jl")
include("phase_retrieval.jl")
include("quadrature.jl")

end
