module Interferometers

using LinearAlgebra

# Simulation functions
export michelson
export sagnac
export fringes
export fringes_gaussian
export gaussian

# Phase recovery functions
export phase_atan
export phase_highgain
export make_cos_first
export quadrature_fit
export quadrature_set
export atan_phase_offset

# filter types multiple dispatch
abstract type specFilter end
struct Gaussian <: specFilter end
struct Square <: specFilter end
export Gaussian, Square


include("simulation_models.jl")
include("phase_retrieval.jl")
include("quadrature.jl")

end
