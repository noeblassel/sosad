using SparseArrays,Arpack,LinearAlgebra,Distributions,ProgressMeter,LaTeXStrings,Plots
include("eigen.jl")
include("pot_opt.jl")


ab_pairs = [(0.5,0.3),(1.0,-0.3),(0.0,0.0)]
markershapes=[:circle,:xcross,:cross]
styles=[:dashdot,:dot,:dash]
colors = [:green,:blue,:red]
βrange = range(0,16,33)

l2_plot = plot(ylabel="Eigenvalues",xlabel=L"\beta",legend=:bottomright)
V_plot = plot(X,V.(X),xlabel=L"x",ylabel=L"V",label="",color=:black)
err_plot = plot(xlabel=L"\beta",ylabel="Error")

nev = 3
m = -Inf
μ0s,us0_harm = harm_alpha(20;N=N,nev=nev)
λs_min = ν0*(μ0s .- 0.5)

for (i,t)=enumerate(ab_pairs)
    a,b = t
    λs = Vector{Float64}[]
    vline!(V_plot,[X[z0]-a/sqrt(10),X[z1]+b/sqrt(10)],linestyle=styles[i],color=colors[i],label="")
    @showprogress for β = βrange
        Ωind = Ω_f(β,a,b)
        try
            evals, us = eigenproblem(β,N,X,Ωind,V;nev=nev,witten=true)
            push!(λs,evals[1:end])
        catch
            push!(λs,fill(NaN,nev))
        end
    end
    λs = reduce(hcat,λs)

    μ1s,us1_harm= harm_alpha(a*sqrt(ν1/2);N=N,nev=nev)
    μ2s,us2_harm = harm_alpha(b*sqrt(ν2/2);N=N,nev=nev)
    
    λs_harm = vcat(ν1*(μ1s .+ 0.5),ν2*(μ2s .+ 0.5))
    p = sortperm(λs_harm)
    λs_harm = λs_harm[p]
    ix = vcat(fill(1,nev),fill(2,nev),fill(3,nev))[p] # 1=i0, 2=z1, 3=z2
    νs = [ν0,ν1,ν2][ix[p]]

    # println(νs)
    # println(ix)
    # println(λs_harm)
    
    λs_tot = sort(vcat(λs_harm,λs_min))
    global m = max(m,λs_tot[nev])

    j = maximum(i for i=1:nev if λs_harm[i] <= λs_tot[nev])
    λs[nev,λs[1,:] .> 5] .= NaN # omit top eigenvalue if bottom eigenvalue failed to converge
    scatter!(l2_plot,βrange,λs',label="",markershape=markershapes[i],color=colors[i],markerstrokecolor=colors[i],markercolor=((colors[i] == :green) ? :white : colors[i]))
    hline!(l2_plot,λs_harm[1:j],linestyle=styles[i],color=colors[i],label="")

    diffs = abs.(λs_tot[2:j] .- λs[2:j,:])
    println(diffs)
    println(size(diffs)," ",size(βrange))
    scatter!(err_plot,βrange,diffs',color=colors[i],markerstrokecolor=colors[i],markercolor=((colors[i] == :green) ? :white : colors[i]),label="")

end

hline!(l2_plot,λs_min[1:nev],linestyle=:dashdotdot,color=:black,label="",title="")

ylims!(l2_plot,0.0,1.1m)
savefig(l2_plot,"plots/l2.pdf")
title!(V_plot,L"\beta=10")
savefig(V_plot,"plots/V.pdf")
savefig(err_plot,"errors.pdf")