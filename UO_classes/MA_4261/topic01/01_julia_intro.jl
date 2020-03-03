using BenchmarkTools
function mygemm!(C, A, B)
# gemm General Matrix Multiplication
# syntax sugaring !
# C := C + A * B
    M, N = size(C)
    K = size(A,2)

    @assert size(C) == (M,N)
    @assert size(A) == (M,K)
    @assert size(B) == (K,N)


    @inbounds for i = 1:M
        @inbounds for j = 1:N
            @inbounds for p = 1:K
                C[i, j] = C[i, j] + A[i,p] * B[p,j]
            end
        end
    end
    return C
end
# 
# C = rand(10,12)
#
# A = rand(10,15)
#
# B = rand(15,12)
#
# mygemm!(C,A,B)
#
# D = C + A*B
#
# extrema(D - C) # get minimum and maximum value


using LinearAlgebra
let # inside let are global scale variables
    M,N,K = 256, 256, 256
    C = rand(M,N)
    A = rand(M,K)
    B = rand(K,N)

    #@time D = C + A * B
    @time mul!(C,A,B)
    @time mygemm!(C,A,B)

    nsamples = 10

    t_mul! = floatmax(Float64)
    for s = 1:nsamples
        t = @elapsed mul!(C,A,B)
        t_mul! = min(t_mul!, t)
    end
    @show t_mul!

    t_mygemm! = floatmax(Float64)
    for s = 1:nsamples
        t = @elapsed mygemm!(C,A,B)
        t_mygemm! = min(t_mygemm!, t)
    end
    @show t_mygemm!

end
nothing
