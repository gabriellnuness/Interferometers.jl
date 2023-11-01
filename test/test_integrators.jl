using Test
using Interferometers
using PyPlot


@testset "Integrators: dy/dx= sin(2π*x)" begin 
    function fun(y, p, x)
        return sin(2π*x)
    end

    dx = 0.1
    x = 0:dx:1
    n = length(x)
    y_analytical =@. sin(π*x)^2/π
    y_rk4 = zeros(n)
    y_bs3 = zeros(n)
    y_euler = zeros(n)
    y0 = 0
    y_rk4[1] = y0
    y_euler[1] = y0
    y_bs3[1] = y0
    p = 0
    for i =1:n-1
        y_rk4[i+1] = Interferometers.integrate(du=fun, u=y_rk4[i], dt=dx, p=p, t=x[i])
        y_euler[i+1] = Interferometers.integrate(Euler, du=fun, u=y_euler[i], dt=dx, p=p, t=x[i])
        y_bs3[i+1] = Interferometers.integrate(BS3, du=fun, u=y_euler[i], dt=dx, p=p, t=x[i])
    end
    figure()
    plot(x, y_analytical, label="analytical",linewidth=3, alpha=.3,color="black")
    plot(x,y_rk4, label="rk4")
    plot(x,y_bs3, label="bs3")
    plot(x,y_euler, label="euler")
    legend()
end



@testset "Integrators: dy/dx= y*exp(-x)" begin 
    function fun(x, params, t)
        return x*exp(-t)
    end
    dx = 0.6
    x = 0:dx:5
    n = length(x)
    y_analytical =@. exp(1)*exp(-exp(-x))
    y_rk4 = zeros(n)
    y_bs3 = zeros(n)
    y_euler = zeros(n)
    y0 = 1
    y_rk4[1] = y0
    y_euler[1] = y0
    y_bs3[1] = y0
    p = 0
    for i = 1:n-1
        y_euler[i+1] = Interferometers.integrate(Euler, du=fun, u=y_euler[i], dt=dx, p=p, t=x[i])
        y_bs3[i+1] = Interferometers.integrate(BS3, du=fun, u=y_euler[i], dt=dx, p=p, t=x[i])
        y_rk4[i+1] = Interferometers.integrate(du=fun, u=y_rk4[i], dt=dx, p=p, t=x[i])
    end
    figure()
    plot(x, y_analytical, label="analytical",linewidth=3, alpha=.3,color="black")
    plot(x, y_rk4, label="rk4")
    plot(x,y_bs3, label="bs3")
    plot(x,y_euler, label="euler")
    legend()
end






@testset "Integrators: Lorenz attractor" begin 
    function fun(u, params, t)
        x, y, z = u
        a, b, c = params

        dx = a*(y-x)
        dy = x*(c-z) - y
        dz = x*y - b*z

        du = [dx, dy, dz]
        return du
    end
    dt = 0.001
    t = 0:dt:50
    n = length(t)
    u_rk4 = zeros(n,3)
    u_bs3 = zeros(n,3)
    u_euler = zeros(n,3)
    u0 = [1,0,0]
    u_rk4[1,:] = u0
    u_euler[1,:] = u0
    u_bs3[1,:] = u0
    p = (10, 8/3, 28)
    for i = 1:n-1
        u_rk4[i+1,:] = Interferometers.integrate(du=fun, u=u_rk4[i,:], dt=dt, p=p, t=t[i])
        u_euler[i+1,:] = Interferometers.integrate(Euler, du=fun, u=u_euler[i,:], dt=dt, p=p, t=t[i])
        u_bs3[i+1,:] = Interferometers.integrate(BS3, du=fun, u=u_euler[i,:], dt=dt, p=p, t=t[i])
    end

    figure();plot3D(u_rk4[:,1],u_rk4[:,2],u_rk4[:,3],linewidth=.5,label="rk4");legend()
    figure();plot3D(u_bs3[:,1],u_bs3[:,2],u_bs3[:,3],linewidth=.5,label="bs3");legend()
    figure();plot3D(u_euler[:,1],u_euler[:,2],u_euler[:,3],linewidth=.5,label="euler");legend()
end




@testset "Integrators: dy/dx= k(T_out - y)" begin 
    function fun(y, p, x)
        temp_coef, temp_out = p
        return temp_coef*(temp_out - y)
    end
    dx = 0.05
    x = 0:dx:1
    n = length(x)
    y_rk4 = zeros(n)
    y_bs3 = zeros(n)
    y_euler = zeros(n)
    y0 = 100
    y_rk4[1] = y0
    y_euler[1] = y0
    y_bs3[1] = y0
    temp_out = 5
    temp_coef = 10
    p = (temp_coef, temp_out)
    for i =1:n-1
        y_rk4[i+1] = Interferometers.integrate(du=fun, u=y_rk4[i], dt=dx, p=p, t=x[i])
        y_euler[i+1] = Interferometers.integrate(Euler, du=fun, u=y_euler[i], dt=dx, p=p, t=x[i])
        y_bs3[i+1] = Interferometers.integrate(BS3, du=fun, u=y_euler[i], dt=dx, p=p, t=x[i])
    end
    y_analytical =@. temp_out+(y0-temp_out)*exp(-temp_coef*x)
    figure()
    plot(x, y_analytical, label="analytical",linewidth=3, alpha=.3,color="black")
    plot(x,y_rk4, label="rk4")
    plot(x,y_bs3, label="bs3")
    plot(x,y_euler, label="euler")
    legend()

end