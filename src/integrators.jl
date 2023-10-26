function forward_euler(dy, y, dt, p, t=0)
    y_next = y + (dy(y,p,t) * dt)
    return y_next
end

function rk4(dy, y, dt, p, t=0)
    k1 = dy(y,            p,  t)
    k2 = dy(y + dt/2*k1,  p,  t + dt/2)
    k3 = dy(y + dt/2*k2,  p,  t + dt/2)
    k4 = dy(y + dt*k3,    p,  t + dt)
    y_next = y + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    return y_next
end


