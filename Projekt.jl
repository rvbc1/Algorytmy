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


#x = [1; 1]

array_size = 100

A = rand(array_size,array_size)
x = zeros(array_size,1)
b = rand(array_size,1)
#println(b)

#A = [45 23 23; 12 86 2; 14 53 45]
#x = [1; 1; 1]
#b = [-5646; 54355; 43]
#b = [0.14745; 0.406418; 0.624058]



#perf_gres = performance_gmres(A, x, b, 0.1, 1000)
#print_performance(perf_gres)

#println()

perf_bicgstab = performance_bicgstab(A, x, b, 0.1, 100000)
print_performance(perf_bicgstab)
