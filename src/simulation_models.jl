"""
    Low coherence spatial fringes simulation


# parameters:
*   `λ`:                Center wavelength.
*   `Δλ`:               Wavelength bandwidth.
*   `ΔL`:               Optical path lentgh to analyze.
*   `ΔL_res`:           Optional resolution for the OPL.


# returns:
*   `length`:           Physical length.
*   `Interferogram`:    Normalized interferogram intensity.
"""
function fringes(λ, Δλ, ΔL, ΔL_res=1e-6)

    opl_range = range(start=-ΔL, stop=ΔL, step=ΔL_res) 
    wavelength_range = range(start=(λ-Δλ/2), stop=(λ+Δλ/2), length=size(opl_range,1)) 
    interferogram = zero(opl_range)

    for (i, dl) in enumerate(opl_range)
        ϕ =@. 2π * dl / wavelength_range # each λ interference with each distance
        fringe = michelson.(ϕ)
        interferogram[i] = sum(fringe)
    end
    interferogram = interferogram / maximum(interferogram)

    return(length=opl_range, interferogram = interferogram)
end



""" Spectrum input Gaussian filtered """
function fringes_gaussian(λ, Δλ, ΔL, ΔL_res=1e-8)

    Δλ = 2*Δλ
    opl_range = range(start=-ΔL, stop=ΔL, step=ΔL_res) 
    wavelength_range = range(start=(λ-Δλ/2), stop=(λ+Δλ/2), length=size(opl_range,1)) 
    
    gaus_filter = gaussian(wavelength_range, Δλ/2)

    interferogram = zero(opl_range)

    for (i, dl) in enumerate(opl_range)
        ϕ =@. 2π * dl / wavelength_range # each λ interference with each distance
        fringe =@. michelson(ϕ) * gaus_filter
        interferogram[i] = sum(fringe)
    end
    interferogram = interferogram / maximum(interferogram)

    return(length=opl_range, interferogram = interferogram)
end




"""
    Simulation of a Michelson interferometer
"""
function michelson(ϕ,A=1,V=1)
    return A + A*V*cos(ϕ)
end



"""
    Simulation of a Sagnac interferometer
"""
function sagnac(t, τ, ϕs, ϕm)
    0.5 + 0.5*cos(ϕs(t) + ϕm(t) - ϕm(t-τ))
end


"""
    Gaussian normalized function to filter spectrum
"""
function gaussian(x, sigma)

    μ = (x[end]-x[1])/2
    gaus =@. 1/√(2π)/sigma * exp(-1/2 * (x-μ)^2 / (sigma^2))
    gaus = gaus / maximum(gaus)

    return gaus
end