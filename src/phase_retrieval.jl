"""
(phase, spurious_phase) = phase_highgain(arr_cos, arr_sin, dt, gain, e)

This function demodulates the interferometric phase by
using the nonlinear control technique based on sliding-modes.

ref:  FELÃO, Luiz Henrique Vitti. "High gain approach and sliding mode
      control applied to quadrature interferometer". 2019.
      https://repositorio.unesp.br/handle/11449/190782
"""
function phase_highgain(arr_cos::Vector, arr_sin::Vector, dt, gain, e)
    # TODO: review this function and implement tests
    control_phase = zeros(length(arr_cos))
    x1 = zeros(length(arr_cos)-1)
    f = zeros(length(arr_cos)-1)
    u = zeros(length(arr_cos)-1)
    
    for i = 1:length(arr_cos)-1
        
        x1[i]   = arr_cos[i]*cos(control_phase[i]) - arr_sin[i]*sin(control_phase[i])
        f[i]    = arr_cos[i]*sin(control_phase[i]) + arr_sin[i]*cos(control_phase[i])
        u[i]    = -gain * sigmoid(x1[i]*f[i], e)

        if i == 1
           control_phase[i+1] = 0 
        else
           control_phase[i+1] = control_phase[i] + ((u[i-1]+u[i])*dt/2) # Euler integration
        end
            
    end

    phase = -control_phase
    spurious_phase = sum(phase)/length(phase)

    (phase, spurious_phase)

end








"""
sigmoid(arr, e) = @. arr / (abs(arr) + e)

Sigmoid function of a signal.

input: arr

"""
sigmoid(arr, e) = @. arr / (abs(arr) + e)





"""
[output_phase, spurious_phase] = phase_atan(arr_cos, arr_sin)

This funcion demodulates the interferometer phase:
Δϕ = tan⁻¹(sin/cos) + m⋅π  

Inputs:   arr_cos, -> cos     -> Interferometric signals in quadrature without offset
          arr_sin  -> sin 
                       

Outputs:  output_phase        -> Recovered phase
          spurious_phase      -> Spurious phase (fixed value)

ref:  LEMES, "novas configurações de interferômetros de quadrature e de técnicas
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


    # TODO: review theory about offset removal 
    # Before recovering the phase, the algorithm works better
    # works up to Δϕ = @. π*sin(2π*freq_mod*t) + 1.9*π
    # at the bottom it works up to Δϕ = @. π*sin(2π*freq_mod*t) + 0.9*π
    # However it never works with negative phase.

    phase_offset = (maximum(phase) + minimum(phase)) / 2
    turns_phase_offset = phase_offset/(2π)
    
    incomplete_offset_turns = turns_phase_offset - floor(turns_phase_offset)
    
    if incomplete_offset_turns > 0.5
        incomplete_offset_turns = 1-incomplete_offset_turns
        spurious_phase = -incomplete_offset_turns*2π
    else
        spurious_phase = incomplete_offset_turns*2π
    end

    # removing offset level
    final_phase_offset = phase_offset - spurious_phase
    phase = phase .- final_phase_offset
    

    (phase, spurious_phase)
end