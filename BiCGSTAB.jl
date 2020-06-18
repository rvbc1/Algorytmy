using LinearAlgebra

function bicgstab(A, x, b, tolerance, max_itration = 1000)
    r = b - A*x 
    r0 = p = copy(r) 
    rho = dot(r',r0)
    iterations = max_itration;
    for iteration_counter = 1:max_itration
        v = A*p
        alpha = rho/dot(v',r0)
        s = r - alpha*v
        t = A*s
        omega = dot(t',s) / dot(t',t)
        x = x + alpha*p + omega*s
        r_next = s - omega*t

        if(norm(s,2)<tolerance )
            iterations = iteration_counter
            break
        end  

        if(rho==0)
            iterations = iteration_counter
            break
        end

        rho_next = dot(r_next',r0)
        beta = (rho_next / rho) * (alpha/omega)
        p = r_next + beta*(p - omega*v)
        r = r_next
        rho = rho_next
    end
    return x, iterations
end