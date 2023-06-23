module Interferometers

using LinearAlgebra

# Simulation functions
export michelson
export sagnac
# Phase recovery functions
export phase_atan
export phase_highgain
export check_ellipse_rotation
export quadrature_fit
export quadrature_set


include("michelson.jl")
include("sagnac.jl")
include("phase_highgain.jl")
include("phase_atan.jl")
include("quadrature.jl")

end
