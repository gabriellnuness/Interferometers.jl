# Helper functions to analyze data
function fit_line(x,y)
    A = [x ones(length(x))] 
    u = (A'*A)\A'*y

    y_fit = @. u[1]*x + u[2]

    return (y_fit, u)
end
function fit_exp(x,y)
    # y = a⋅eᵇˣ
    y = log.(y)
    sol = fit_line(x,y)
    
    a = exp(sol[2][2])
    b = sol[2][1]
    
    u = [a; b]
    y_fit = @. a*exp(x*b)
    return (y_fit, u)
end
function fit_power(x,y)
    # y = a⋅xᵇ
    x1 = log.(x)
    y1 = log.(y)
    sol = fit_line(x1,y1)
    
    a = exp(sol[2][2])
    b = sol[2][1]
    u = [a; b]
    y_fit = a.*x.^b
    return (y_fit, u)
end