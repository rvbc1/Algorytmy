using LinearAlgebra

function gmres(A, x, b, tolerance, max_iterations)
    m,n = size(A,1), size(A,2)
    r = b - A*x

    b_norm = norm(b)
    error = norm(r) / b_norm

    sn = zeros(max_iterations, 1)
    cs = zeros(max_iterations, 1)
    e1 = zeros(max_iterations, 1)
    e1[1] = 1
    e = [error]
    r_norm = norm(r)
    Q = zeros(n, max_iterations)
    Q[:,1] = r / r_norm
    beta = r_norm * e1

    H = zeros(1, 0)
    k_end = 1

    for k = 1:max_iterations
        #reshape matrix H
        H_add = zeros(k, 1)
        H = hcat(H,H_add)
        H_add = zeros(1, k)
        H = vcat(H,H_add)

        (H[1:k+1, k], Q[:, k+1]) = arnoldi(A, Q, k)
        (H[1:k+1, k], cs[k], sn[k]) = apply_givens_rotation(H[1:k+1,k], cs, sn, k)

        beta[k + 1] = -sn[k] * beta[k]
        beta[k] = cs[k] * beta[k]
        error = abs(beta[k + 1]) / b_norm

        e = [e; error]

        if (error <= tolerance)
            print("in tolerance afer ")
            print(k)
            println(" iterations")
            k_end = k
            break
        end
    end
    y = H[1:k_end, 1:k_end] \ beta[1:k_end]
    x = x + Q[:, 1:k_end] * y

    return x
end

function arnoldi(A, Q, k)
    q = A*Q[:,k]
    h = zeros(k + 1)
    for i = 1:k
        h[i] = q' * Q[:, i]
        q = q - h[i] * Q[:, i]
    end
    h[k + 1] = norm(q)
    q = q / h[k + 1]
    return (h, q)
end

function apply_givens_rotation(h, cs, sn, k)
    for i = 1:k-1
        temp   =  cs[i] * h[i] + sn[i] * h[i + 1]
        h[i+1] = -sn[i] * h[i] + cs[i] * h[i + 1]
        h[i]   = temp
    end
  
    (cs_k, sn_k) = givens_rotation(h[k], h[k + 1])
  
    h[k] = cs_k * h[k] + sn_k * h[k + 1]
    h[k + 1] = 0.0
    return (h, cs_k, sn_k)
end

function givens_rotation(v1, v2)
    if (v1 == 0)
        cs = 0
        sn = 1
    else
        t = sqrt(v1^2 + v2^2)
        cs = abs(v1) / t
        sn = cs * v2 / v1
    end
    return (cs, sn)
end

A = [1 1 -3 1; -5 3 -4 1; 1 0 2 -1; 1 2 0 0]
x = [1; 1; 1; 1]
b = [2; 0; 1; 12]
x0 = gmres(A, x, b, 0.1, 1000)

println(x0)
println(A*x0)