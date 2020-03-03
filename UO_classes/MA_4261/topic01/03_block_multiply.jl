include("gemm.jl")

using LinearAlgebra

let
    nsamples = 10
    ns = 48:48:1500
    for iter = 1:length(ns)
        n = m = k = ns[iter]
        C = rand(n,m)
        A = rand(n,k)
        B = rand(k,m)

        # julia mul!
        t_mul! = floatmax(Float64)
        for s = 1:nsamples
            t = @elapsed mul!(C,A,B)
            t_mul! = min(t,t_mul!)
        end
        GFLOPS_mul![iter] = 2*n*m*k / t_mul! / 1e9

        t_ijp! = floatmax(Float64)
        for s = 1:nsamples
            t = @elapsed mygemm_ijp!(C,A,B)
            t_ijp! = min(t,t_ijp!)
        end

    end
    df = DataFrame(m = ns,

    )
end
