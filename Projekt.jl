include("BiCGSTAB.jl") 
include("GMRES.jl")

using Printf

function performance_bicgstab(A, x0, b, tolerance, max_iterations)
    timed = @timed x = bicgstab(A, x0, b, tolerance, max_iterations)
    exe_time = timed[2]
    used_mem = timed[3]
    iterations = x[2]

    return exe_time, used_mem, iterations
end

function performance_gmres(A, x0, b, tolerance, max_iterations)
    timed = @timed x = gmres(A, x0, b, tolerance, max_iterations)
    exe_time = timed[2]
    used_mem = timed[3]
    iterations = x[2]

    return exe_time, used_mem, iterations
end

function print_performance(performance)
    println("Execution performance:")
    print("\t Time: ")
    print(performance[1])
    println("s")
    print("\t Memory: ")   
    if (performance[2] < (1024))
        print(performance[2])
        println(" B")
    elseif (performance[2] < (1024*1024))
        @printf "%.3f" performance[2]/1024
        #print(performance[2]/1024)
        println(" KiB")
    elseif (performance[2] < (1024*1024*1024))
        @printf "%.3f" performance[2]/1024/1024
        #print(performance[2]/1024/1024)
        println(" MiB")
    else
        @printf "%.3f" performance[2]/1024/1024/1024
        #print(performance[2]/1024/1024/1024)
        println(" GiB")
    end
    print("\t Iterations: ")  
    println(performance[3])
end


array_size = 1000
for i = 1:10

    println()
    print(i)
    println(":")
    A = rand(array_size,array_size)
    print("Cond A: ")
    println(cond(A))
    x = zeros(array_size,1)
    b = rand(array_size,1)

    tolerance = 0.1

    max_itration = 10000




    perf_gres = performance_gmres(A, x, b, tolerance, max_itration)
    print_performance(perf_gres)

    println()

    perf_bicgstab = performance_bicgstab(A, x, b, tolerance, max_itration)
    print_performance(perf_bicgstab)
end
