"""
[signal_1_cos, signal_2_sin, phase, gain_ratio, offset_1, offset_2] = fit_quadrature(signal_1, signal_2)

This function fits an ellipse to the Lissajous curve of signal_1 versus signal_2
and returns as output signals in quadrature (90 deg between them,
a Lissajous circle).
*Obs: The precision of this fitting impacts the phase retrieval error.
  
inputs:
    signal_1, signal_2              -> Interferometric signals out of phase [column vector]

outputs:
    signal_1_cos, signal_2_sin      -> Interferometric signals [column vectors] in quadrature (90 deg)
    phase                           -> Phase difference from quadrature (in rad)
    gain_ratio                      -> Gain proportion between signal_1 and signal_2, gain_ratio = A1*V1/(A2*V2);
    offset_1, offset_2              -> Offsets of signal_1 and signal_2.

ref: "HEYDEMANN, P. L. M. Determination and correction of 
quadrature fringe measurement errors in interferometers. Applied Optics. 
New York, v. 20, n 19, p. 3382-3384, 1981."
"""
function quadrature_fit(signal_1::Vector, signal_2::Vector)
      
    # thesis: eq (49)
    # A*signal_1^2 + B*signal_2^2 + C*signal_1*signal_2 + D*signal_1 + E*signal_2 = 1
    @assert size(signal_1) == size(signal_2)

    X = hcat(signal_1.^2, signal_2.^2, signal_1.*signal_2, signal_1, signal_2)

    I = ones(length(signal_1)) #
    
    # b = [A; B; C; D; E]
    b = (X'*X)\(X')*I         

    # Thesis: eq (61)
    A = b[1]; B = b[2]; C = b[3]; D = b[4]; E = b[5] 
   
    gain_ratio     = sqrt(B/A)
    phase = asin(C/(2*gain_ratio*A))
    offset_1     = (2*B*D-E*C) / (C^2-4*A*B)
    offset_2     = (2*A*E-D*C) / (C^2-4*A*B)
    
    (phase, gain_ratio, offset_1, offset_2)
end


"""
[cos_signal sin_signal] = quadrature_set(alpha, r, p, q)

This function is a subset of fit_quadrature.m function
the inputs are the fitted values necessary to correct the
quadrature.

Input: 
    phase                           -> Phase difference from quadrature (in rad)
    gain_ratio                      -> Gain proportion between signal_1 and signal_2, gain_ratio = A1*V1/(A2*V2);
    offset_1, offset_2              -> Offsets of signal_1 and signal_2.

Output: cosine signal,
            sine signal.

ref: Felão, "High gain approach and sliding mode control
      applied to quadrature interferometer", Thesis, UNESP,
      2019. https://repositorio.unesp.br/handle/11449/190782.
"""
function quadrature_set(signal_1, signal_2, phase, gain_ratio, offset_1, offset_2)

    # Thesis: eq. (43) and (44)
    cos_signal = @. signal_1 - offset_1;
    sin_signal = @. (cos_signal*sin(phase) + gain_ratio*(signal_2-offset_2)) / cos(phase);

    (cos_signal, sin_signal)
end


"""
(signal_1_cos, signal_2_sin) = check_ellipse_rotation(signal_1, signal_2, δ)

Check if signal is representative of a cosine or sine
Make the first signal the cosine version.
The calculation is done through the rotation direction of the Lissajous figure
if x: cos and y: sin, than the direction is counter clockwise.
The rotation direction is measured by the cross product of vectors.

"""
function check_ellipse_rotation(x::Vector, y::Vector, δ::Int)
    
    p1 = 1
    p2 = p1 + δ
    p3 = p2 + δ
    vector_1 = [x[p1]-x[p2]; y[p1]-y[p2]; 0]
    vector_2 = [x[p2]-x[p3]; y[p2]-y[p3]; 0]
    vec_prod = cross(vector_1, vector_2)

    if vec_prod[3] <= 0
    # z-vector entering the screen
        x_cos = y
        y_sin = x

        println("Inverted by check_ellipse_rotation()")
        return (x_cos, y_sin)
    end
    
    println("Not inverted by check_ellipse_rotation()")
    (x, y)
end