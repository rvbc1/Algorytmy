using LinearAlgebra

function bicgstab()
    #A = 4;
    A = [1 2; 3 4]
    x = [1; 1]
    #x = 6;
    #b = 3;
    b = [-39; 30]
    max_itration = 10;

    normal_to_b = norm( b );

    r = b - A*x 
    r0 = p = r 
    

    for iter = 1:max_itration
        alpha = (r'*r0)/((A*p)'*r0)
        s = r - alpha*A*p
        omega = ((A*s)'*s) / ((A*s)'*(A*s))
        x = x + alpha*p + omega*s
        r_next = s - omega*(A*s)

        error = norm(r) / normal_to_b; 
        if(error < 0.1)
            print("in tolerance afer ")
            print(iter)
            println(" iterations")
            break
        end

        if ( omega == 0.0 )
            println("zero omega")
            break
        end

        beta = ((r_next'*r0) / (r'*r0)) * (alpha/omega)
        p = r_next + beta*(p - omega*(A*p))
        r = r_next
    end
    return x
end

x0 = bicgstab()

A = [1 2; 3 4]


println(x0)
println(A*x0)