#!/bin/env julia
#OAR -l /nodes=1/core=8,walltime=25

using SparseArrays,Arpack,LinearAlgebra,Distributions,ProgressMeter,DelimitedFiles
include("eigen.jl")
include("pot_opt.jl")

# β = parse(Int,ARGS[1])
β = 10

function main(β::Int,N_bins=200)
    minf = -1.0
    inf = 1.5
    Φ(x) = cdf(Normal(),Float64(x))
    a_range = range(minf,inf,N_bins)
    b_range = range(minf,inf,N_bins)

    λ1s_sc = fill(NaN,(N_bins,N_bins))
    λ2s_sc = fill(NaN,(N_bins,N_bins))

    μ1_sc = fill(NaN,N_bins)
    μ2_sc = fill(NaN,N_bins)

    if β == -1 # compute semi-classical quantities (fast)
        @showprogress for i = eachindex(a_range)
            a = a_range[i]
            μ1s,w1s = harm_alpha(a*sqrt(ν1/2);N=N,nev=1)
            μ1_sc[i] = first(μ1s)
        end

        @showprogress for j = eachindex(b_range)
            b = b_range[j]
            μ2s,w2s = harm_alpha(b*sqrt(ν2/2);N=N,nev=1)
            μ2_sc[j] = first(μ2s)
        end

        writedlm("m1_harm.csv",μ1_sc)
        writedlm("m2_harm.csv",μ2_sc)

        @showprogress for i = eachindex(a_range)
            for j = eachindex(b_range)
                a = a_range[i]
                b = b_range[j]

                μ1 = μ1_sc[i]
                μ2 = μ2_sc[j]

                λ2 = minimum([ν1*(μ1+0.5),ν2*(μ2+0.5),ν0])
                λ1 = (sqrt(ν1)/Φ(sqrt(ν1)*a) + sqrt(ν2)/Φ(sqrt(ν2)*b))*sqrt(ν0)/2π

                λ1s_sc[i,j] = λ1
                λ2s_sc[i,j] = λ2

            end
        end

        writedlm("l1_sc.csv",λ1s_sc)
        writedlm("l2_sc.csv",λ2s_sc)

    else # compute finite-dimensional quantities (slow)

        λ1_refs = Float64[]
        λ2_refs = Float64[]

        λ1s = fill(NaN,(N_bins,N_bins))
        λ2s = fill(NaN,(N_bins,N_bins))

        Ωind_ = Ω_f(β,0,0)
        λs_,us_,L = eigenproblem(β,N,X,Ωind_,V;nev=2,witten=true)
        λ1_,λ2_ = λs_

        if(λ1_ > 5)
            λ2_ = λ1_
            λ1_ = NaN
        end

        l1_ref = λ1_
        l2_ref = λ2_

        @showprogress for i = eachindex(a_range)
            for j = eachindex(b_range)
                a = a_range[i]
                b = b_range[j]

                Ωind = Ω_f(β,a,b)

                if length(Ωind) > 100
                    λs, us = eigenproblem(β,N,X,Ωind,V;nev=2,witten=true,L=L)

                    λ1,λ2 = λs

                    if(λ1 > 5)
                        λ2 = λ1
                        λ1 = NaN
                    end
                    
                    λ1s[i,j] = λ1 / l1_ref
                    λ2s[i,j] = λ2 / l2_ref

                end
            end
        end

        writedlm("l1_red_$(β).csv",λ1s)
        writedlm("l2_red_$(β).csv",λ2s)
    end

end

# main(β)