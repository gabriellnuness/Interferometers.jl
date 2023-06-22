"""
[cos_signal sin_signal] = quadrature_set(alpha, r, p, q)

This function is a subset of fit_quadrature.m function
the inputs are the fitted values necessary to correct the
quadrature.

Input:    alpha,
              r,
              p,
              q -> ellipse fitting constant values [1,1]

Output: cosine signal,
          sine signal.

ref: Felão, "High gain approach and sliding mode control
      applied to quadrature interferometer", Thesis, UNESP,
      2019. https://repositorio.unesp.br/handle/11449/190782.
"""
function quadrature_set(signal_1, signal_2, alpha, r, p, q)

    # Thesis: eq. (43) and (44)
    cos_signal = @. signal_1 - p;
    sin_signal = @. (cos_signal*sin(alpha) + r*(signal_2-q)) / cos(alpha);

    (cos_signal, sin_signal)
end