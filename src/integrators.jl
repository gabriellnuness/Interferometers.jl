function integrate(::Type{Euler}; du, u, dt, p=nothing, t=0)
    u_next = u + (du(u,p,t) * dt)
    return u_next
end

function integrate(::Type{BS3}=BS3; du, u, dt, p=nothing, t=0)
    k1 = du(u,p,t)
    k2 = du(u+1/2*dt*k1, p, t+1/2*dt)
    k3 = du(u+3/4*dt*k2, p, t+3/4*dt)
    u_next = u + 2/9*dt*k1 + 1/3*dt*k2 + 4/9*dt*k3
    k4 = du(u_next, p, t+dt)
    z_next = u + dt*(7/24*k1 + 1/4*k2 + 1/3*k3 +1/8*k4)
    return z_next
end

function integrate(::Type{RK4}; du, u, dt, p=nothing, t=0)
    k1 = du(u,p,t)
    k2 = du(u+dt/2*k1, p, t+dt/2)
    k3 = du(u+dt/2*k2, p, t+dt/2)
    k4 = du(u+dt*k3, p, t+dt)
    u_next = u + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    return u_next
end


