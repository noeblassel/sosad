using Arpack,SparseArrays,DelimitedFiles,ProgressMeter,Plots

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


function main(nev=4,nbins=20)
    ξ1min = -0.5

    imin,imax = 1,132
    jmin,jmax = 870,999

    bounds = [(i,j) for i=imin:imax, j=jmin:jmax]

    ε = 0.5
    β = 2.0

    L1s = zeros(nbins,nbins)
    L2s = zeros(nbins,nbins)


    M = readdlm("./series/FD_quad_eps=0.5_beta=2.0.out")

    ξ1 = M[:,1]
    F1 = M[:,2]
    a1 = M[:,3]

    Jmax = -1
    istar = -1
    jstar = -1
            
    L = comp_effective_generator(β,ξ1,a1,F1)

    @showprogress for (i,j)=bounds
                λs,us = comp_evs(L,ξ1,ξ1[i],ξ1[j],2)
                if λs[2]/λs[1] > Jmax
                    Jmax = λs[2]/λs[1]
                    istar = i
                    jstar = j
                end

                # println(Jmax)
    end

    println(ξ1[istar]," ",ξ1[jstar]," ",Jmax)

    return nothing
end

main()