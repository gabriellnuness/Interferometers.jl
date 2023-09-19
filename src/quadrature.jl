"""
`(phase, gain_ratio, dc1, dc2) = fit_quadrature(arr1, arr2)`

Fits an ellipse to the Lissajous curve of `arr1` versus `arr2`
and returns as output signals in quadrature (90°  between them, i.e., a circle).
**Obs: The precision of this fitting impacts the phase retrieval error.**
  
# parameters:
*    `arr1`, `arr2`:    Interferometric signals out of phase [column vector]

# returns:
*    `phase`:           Phase difference from quadrature (in rad)
*    `gain_ratio`:      Gain proportion between `arr1` and `arr2`, gain_ratio = A₁V₁/(A₂V₂);
*    `dc1`, `dc2`:      Offsets of `arr1` and `arr2`.

# ref: 
    "HEYDEMANN, P. L. M. Determination and correction of 
    quadrature fringe measurement errors in interferometers. Applied Optics. 
    New York, v. 20, n 19, p. 3382-3384, 1981."
"""
function quadrature_fit(arr1::Vector, arr2::Vector)
      
    # thesis: eq (49)
    # A*arr1^2 + B*arr2^2 + C*arr1*arr2 + D*arr1 + E*arr2 = 1
    @assert size(arr1) == size(arr2)

    X = hcat(arr1.^2, arr2.^2, arr1.*arr2, arr1, arr2)

    I = ones(length(arr1)) #
    
    # b = [A; B; C; D; E]
    b = (X'*X)\(X')*I         

    # Thesis: eq (61)
    A = b[1]; B = b[2]; C = b[3]; D = b[4]; E = b[5] 
   
    gain_ratio = sqrt(B/A)
    phase = asin(C/(2*gain_ratio*A))
    dc1     = (2*B*D-E*C) / (C^2-4*A*B)
    dc2     = (2*A*E-D*C) / (C^2-4*A*B)

    (phase=phase, gain_ratio=gain_ratio, offset1=dc1, offset2=dc2)
end


"""
`quadrature_set(arr1, arr2, phase, gain_ratio, dc1, dc2)`

Modifies the signals to be in quadrature based on the values obtained previously from the fit.

# parameters: ellipse fitted values
*   `phase`: Phase difference from quadrature in rad. 
*   `gain_ratio`: Gain proportion between signal_1 and signal_2, gain_ratio = A₁V₁/(A₂V₂);
*   `dc1`, `dc2`: Offsets of signal_1 and signal_2.

# returns: named Tuple with corrected signals
* cosine signal,
* sine signal.

# ref:
    Felão, "High gain approach and sliding mode control
    applied to quadrature interferometer", Thesis, UNESP,
    2019. https://repositorio.unesp.br/handle/11449/190782.
"""
function quadrature_set(arr1, arr2, phase, gain_ratio, dc1, dc2)

    # Thesis: eq. (43) and (44)
    cos_signal = @. arr1 - dc1;
    sin_signal = @. 1/cos(phase)*(cos_signal*sin(phase) + gain_ratio*(arr2-dc2));

    (arr_cos = cos_signal, arr_sin = sin_signal)
end
function quadrature_set(arr1, arr2, quad)
    return quadrature_set(arr1, arr2, quad[1], quad[2], quad[3], quad[4])
end


"""
`make_cos_first(arr1, arr2, δ)`

Check if signal is representative of a cosine or sine
Make the first signal the cosine version.

# parameters:
*   `arr1`, `arr2`: input arrays.
*   `δ`: distance between points of the Lissajous curve to take the vectors.

The calculation is done through the rotation direction of the Lissajous figure
if x: cos and y: sin, than the rotation direction is counter clockwise.
The rotation direction is measured by the cross product of vectors.

"""
function make_cos_first(arr1::Vector, arr2::Vector, δ::Int)
    
    p1 = 1
    p2 = p1 + δ
    p3 = p2 + δ
    vector_1 = [arr1[p1]-arr1[p2]; arr2[p1]-arr2[p2]; 0]
    vector_2 = [arr1[p2]-arr1[p3]; arr2[p2]-arr2[p3]; 0]
    vec_prod = cross(vector_1, vector_2)

    if vec_prod[3] <= 0
    # z-vector entering the screen
        arr1_cos = arr2
        arr2_sin = arr1

        println("Inverted by make_cos_first()")
        return (arr_cos = arr1_cos, arr_sin = arr2_sin)
    end
    
    println("Not inverted by make_cos_first()")
    (arr_cos = arr1, arr_sin = arr2)
end