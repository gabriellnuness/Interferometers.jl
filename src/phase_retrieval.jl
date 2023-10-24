"""
`phase_highgain(arr_cos, arr_sin, τ, gain, e)`

Demodulates the interferometric phase by
using the nonlinear control technique based on sliding-modes.

# parameters: 
*   `arr_cos`, `arr_sin`:   interferometric signals in quadrature without offset.
*   `τ`:                    data sample period.
*   `gain`:                 sliding modes control gain.
*   `e`:                    sigmoid coefficient (e=0 -> signal function).                       

# returns: named Tuple (phase, offset)
*   `phase`:       Recovered phase
*   `offset`:      Spurious offset phase (ϕ₀)

# ref: 
    FELÃO, Luiz Henrique Vitti. "High gain approach and sliding mode
    control applied to quadrature interferometer". 2019.
    https://repositorio.unesp.br/handle/11449/190782
"""
function phase_highgain(arr_cos::Vector, arr_sin::Vector, τ, gain, e=0)

    ϕc = zeros(length(arr_cos))
    x1 = zeros(length(arr_cos)-1)
    f = zeros(length(arr_cos)-1)
    u = zeros(length(arr_cos)-1)
    
    for i = 1:length(arr_cos)-1
        
        x1[i] = arr_cos[i]*cos(ϕc[i]) - arr_sin[i]*sin(ϕc[i])
        f[i]  = arr_cos[i]*sin(ϕc[i]) + arr_sin[i]*cos(ϕc[i])
        u[i]  = -gain * sigmoid(x1[i]*f[i], e)

        if i == 1
           ϕc[i+1] = 0   # initial condition
        else
           ϕc[i+1] = ϕc[i] + ((u[i-1]+u[i]) * τ/2) # Trapz
        end
            
    end

    ϕ = -ϕc
    spurious_phase = sum(ϕ)/length(ϕ)

    (phase = ϕ, offset = spurious_phase)

end








"""
sigmoid(x, e) = @. x / (abs(x) + e)

Sigmoid function of a signal.

#parameters: 
* `x`: input value
* `e`: sigmoid approximation
"""
function sigmoid(x, e)
    if x==0
        return 0
    end
    return (x / (abs(x) + e))
end





"""
`phase_atan(arr_cos, arr_sin)`

Demodulates the interferometer phase with the arctangent method
Δϕ = tan⁻¹(sin/cos) + m⋅π  

It is limited to π/2 between two samples instead of π by atan2(), 
check unwrap() from DSP package.

# parameters: 
 *   `arr_cos` and `arr_sin`: Interferometric signals in quadrature without offset.
                       

# returns: named Tuple (phase, offset)
*   `phase`:       Recovered phase
*   `offset`:      Spurious offset phase (ϕ₀)

# ref:  
    LEMES, "novas configurações de interferômetros de quadrature e de técnicas
    de detecção de fase óptica baseadas em phase unwrapping",
    dissertation, UNESP, 2014
    https://repositorio.unesp.br/handle/11449/111112
"""
function phase_atan(arr_cos::Vector, arr_sin::Vector)
        
    @assert size(arr_cos) == size(arr_sin)    

    # Angle recovered with discontinuities
    phase = @. abs(atan(arr_sin/arr_cos))    
  
    # Initializing variables
    previous_quadrant = 0
    trig_circle_turns = 0
    
    # Check quadrant of each point and 
    # add to the number of trigonometric circle turns
    for i = 1:length(arr_sin)
        
        if arr_cos[i] >= 0 && arr_sin[i] >= 0
        # 1st
            current_quadrant = 1
        elseif arr_cos[i] <= 0 && arr_sin[i] >= 0
        # 2nd
            current_quadrant = 2
            phase[i] = π - phase[i]
        elseif arr_cos[i] <= 0 && arr_sin[i] <= 0
        # 3rd
            current_quadrant = 3
            phase[i] = π + phase[i]
        elseif arr_cos[i] >= 0 && arr_sin[i] <= 0
        # 4th
            current_quadrant = 4
            phase[i] = 2π - phase[i]
        end  
        
        # counting number of circle turns
        if current_quadrant == 1 && previous_quadrant == 4
            trig_circle_turns = trig_circle_turns + 1
        elseif current_quadrant == 4 && previous_quadrant == 1
            trig_circle_turns = trig_circle_turns - 1
        end
        previous_quadrant = current_quadrant
    
        # Total phase calculation in degrees
        phase[i] = phase[i] + 2π*trig_circle_turns

    end 
    
    (offset, spurious_phase) = atan_phase_offset(phase)
    phase = phase .- offset
    
    return (phase=phase, offset=spurious_phase)
end



function atan_phase_offset(phase)
        offset = (maximum(phase) + minimum(phase)) / 2
        turns = offset/(2π)      
        incomplete_turns = turns - floor(turns)
        
        if incomplete_turns > 0.5
            incomplete_turns = 1 - incomplete_turns # complement
            spurious_phase = -incomplete_turns * 2π
        else
            spurious_phase = incomplete_turns * 2π
        end

       final_phase_offset = offset - spurious_phase

       return (final_phase_offset, spurious_phase)
    end



"""
Maximum value between two phase samples that can still be unwrapped
ϕ < fs⋅π/(2πf)
"""
function check_arctangent_lim_phase(fs,f;lim2sample=π)
    maximum_phase =  fs*lim2sample/(2π*f)
    return maximum_phase
end
function check_arctangent_lim_freq(ϕ,f;lim2sample=π)
    min_fs =  2π*f*ϕ/lim2sample
    return min_fs
end