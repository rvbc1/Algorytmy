using LinearAlgebra

function bicgstab(A, x, b, tolerance, max_itration = 1000)
    norm_b = norm(b);
    r = b - A*x 
    r0 = p = r 
    rho = r'*r0
    iterations = max_itration;
    for iteration_counter = 1:max_itration
        v = A*p
        alpha = rho/(v'*r0)
        s = r - alpha*v
        t = A*s
        omega = (t'*s) / (t'*t)
        x = x + alpha*p + omega*s
        r_next = s - omega*t

        error = norm(r) / norm_b; 
        if(error <= tolerance)
            #print("in tolerance afer ")
            #print(iteration_counter)
            #println(" iterations")
            iterations = iteration_counter
            break
        end

        if ( omega == 0.0 )
            println("zero omega")
            break
        end

        rho_next = (r_next'*r0)
        beta = (rho_next / rho) * (alpha/omega)
        p = r_next + beta*(p - omega*v)
        r = r_next
        rho = rho_next
    end
    return x, iterations
end



