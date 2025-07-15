using Arpack,SparseArrays,DelimitedFiles,ProgressMeter

"""
    ξ : bin centers
    F : free energy
    a : effective diffusion
"""
function comp_effective_generator(β,ξ,a,F)
    N = length(ξ)
    dξ = (last(ξ)-first(ξ))/(N-1)

    @inline prev(i) = mod(i-1,1:N)
    @inline next(i) = mod(i+1,1:N)
    rate(i,j) = exp(-β*(F[j]-F[i])/2)*(a[i]+a[j])/2
    rows = vcat(([i,i,i] for i=1:N)...)
    cols = vcat(([i,prev(i),next(i)] for i=1:N)...)

    factor = inv(β*dξ^2)
    vals = vcat((factor*[rate(i,prev(i))+rate(i,next(i)),-rate(i,prev(i)),-rate(i,next(i))] for i=1:N)...)

    L = sparse(rows,cols,vals)

    return L
end

"""
   Compute Dirichlet eigenvalues of effective generator
"""
function comp_evs(L,ξ,ξmin,ξmax,nev)
    Ωind = (ξmin .< ξ .< ξmax)
    λs,us =  eigs(L[Ωind,Ωind],sigma=0.0,nev=nev,which=:LR,tol=1e-12,maxiter=10000)

    return real.(λs),real.(us)
end


function main(nev=4,nbins=1000)
    βs = [2.0]
    εs = [2.0,1.0,0.5,0.1]

    ξ1min = -0.5
    b1min, b1max = 0.2,0.8

    ξ2min = -2.1
    b2min, b2max = -0.5,0.5

    for β=βs
        for ε=εs
            λ_ξ1s = Vector{Float64}[]
            b1s = Float64[]

            λ_ξ2s = Vector{Float64}[]
            b2s = Float64[]

            M = readdlm("./series/FD_quad_eps=$(ε)_beta=$(β).out")
            ξ1 = M[:,1]
            F1 = M[:,2]
            a1 = M[:,3]
            dξ1 = 0.002

            ξ2 = M[:,4]
            F2 = M[:,5]
            a2 = M[:,6]
            dξ2 = 0.005
            

            L_ξ1 = comp_effective_generator(β,ξ1,a1,F1)
            L_ξ2 = comp_effective_generator(β,ξ2,a2,F2)

            println(size(L_ξ1))
            println(size(L_ξ2))

            b1 = collect(b1min:dξ1:b1max)
            b2 = collect(b2min:dξ2:b2max)

            println(length(b1))
            println(length(b2))

            @showprogress for b=b1
                evs,_ = comp_evs(L_ξ1,ξ1,ξ1min,b,nev)
                push!(λ_ξ1s,evs)
            end

            @showprogress for b=b2
                evs,_ = comp_evs(L_ξ2,ξ2,ξ2min,b,nev)
                push!(λ_ξ2s,evs)
            end

            writedlm("./series/evs_eff_xi1_eps=$(ε)_beta=$(β)_a1=$(ξ1min).out",[b1 hcat(λ_ξ1s...)'])
            writedlm("./series/evs_eff_xi2_eps=$(ε)_beta=$(β)_a1=$(ξ2min).out",[b2 hcat(λ_ξ2s...)'])

        end
    end
    return nothing
end

main()