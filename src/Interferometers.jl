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

# Integrator for sliding mode multiple dispatch
abstract type IntSolver end
struct RK4 <: IntSolver end
struct BS3 <: IntSolver end
struct Euler <: IntSolver end
export RK4, BS3, Euler


include("simulation_models.jl")
include("phase_retrieval.jl")
include("quadrature.jl")
include("integrators.jl")

end
