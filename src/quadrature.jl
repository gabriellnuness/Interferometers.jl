"""
(phase, gain_ratio, dc_1, dc_2) = fit_quadrature(arr_1, arr_2)

This function fits an ellipse to the Lissajous curve of arr_1 versus arr_2
and returns as output signals in quadrature (90 deg between them,
a Lissajous circle).
*Obs: The precision of this fitting impacts the phase retrieval error.
  
inputs:
    arr_1, arr_2              -> Interferometric signals out of phase [column vector]

outputs:
    phase                     -> Phase difference from quadrature (in rad)
    gain_ratio                -> Gain proportion between arr_1 and arr_2, gain_ratio = A1*V1/(A2*V2);
    dc_1, dc_2                -> Offsets of arr_1 and arr_2.

ref: "HEYDEMANN, P. L. M. Determination and correction of 
quadrature fringe measurement errors in interferometers. Applied Optics. 
New York, v. 20, n 19, p. 3382-3384, 1981."
"""
function quadrature_fit(arr_1::Vector, arr_2::Vector)
      
    # thesis: eq (49)
    # A*arr_1^2 + B*arr_2^2 + C*arr_1*arr_2 + D*arr_1 + E*arr_2 = 1
    @assert size(arr_1) == size(arr_2)

    X = hcat(arr_1.^2, arr_2.^2, arr_1.*arr_2, arr_1, arr_2)

    I = ones(length(arr_1)) #
    
    # b = [A; B; C; D; E]
    b = (X'*X)\(X')*I         

    # Thesis: eq (61)
    A = b[1]; B = b[2]; C = b[3]; D = b[4]; E = b[5] 
   
    gain_ratio     = sqrt(B/A)
    phase = asin(C/(2*gain_ratio*A))
    dc_1     = (2*B*D-E*C) / (C^2-4*A*B)
    dc_2     = (2*A*E-D*C) / (C^2-4*A*B)
    
    (phase, gain_ratio, dc_1, dc_2)
end


"""
[cos_signal sin_signal] = quadrature_set(alpha, r, p, q)

This function is a subset of fit_quadrature.m function
the inputs are the fitted values necessary to correct the
quadrature.

Input: 
    phase           -> Phase difference from quadrature (in rad)
    gain_ratio      -> Gain proportion between signal_1 and signal_2, gain_ratio = A1*V1/(A2*V2);
    dc_1, dc_2      -> Offsets of signal_1 and signal_2.

Output: cosine signal,
            sine signal.

ref: Felão, "High gain approach and sliding mode control
      applied to quadrature interferometer", Thesis, UNESP,
      2019. https://repositorio.unesp.br/handle/11449/190782.
"""
function quadrature_set(arr_1, arr_2, phase, gain_ratio, dc_1, dc_2)

    # Thesis: eq. (43) and (44)
    cos_signal = @. arr_1 - dc_1;
    sin_signal = @. 1/cos(phase)*(cos_signal*sin(phase) + gain_ratio*(arr_2-dc_2));

    (cos_signal, sin_signal)
end


"""
(arr_1_cos, arr_2_sin) = make_cos_first(arr_1, arr_2, δ)

Check if signal is representative of a cosine or sine
Make the first signal the cosine version.

The calculation is done through the rotation direction of the Lissajous figure
if x: cos and y: sin, than the rotation direction is counter clockwise.
The rotation direction is measured by the cross product of vectors.

"""
function make_cos_first(arr_1::Vector, arr_2::Vector, δ::Int)
    
    p1 = 1
    p2 = p1 + δ
    p3 = p2 + δ
    vector_1 = [arr_1[p1]-arr_1[p2]; arr_2[p1]-arr_2[p2]; 0]
    vector_2 = [arr_1[p2]-arr_1[p3]; arr_2[p2]-arr_2[p3]; 0]
    vec_prod = cross(vector_1, vector_2)

    if vec_prod[3] <= 0
    # z-vector entering the screen
        arr_1_cos = arr_2
        arr_2_sin = arr_1

        println("Inverted by make_cos_first()")
        return (arr_1_cos, arr_2_sin)
    end
    
    println("Not inverted by make_cos_first()")
    (arr_1, arr_2)
end